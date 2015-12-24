from doit.cmd_base import ModuleTaskLoader
from doit.compat import is_bound_method
import inspect
from doit import loader
from doit.exceptions import InvalidDodoFile
import six
from doit.loader import generate_tasks, TASK_STRING
from .parser import config_from_args

class WrapitLoader(ModuleTaskLoader):
    """Loader class that does the loading of tasks"""
    def __init__(self, args, mod_dict):
        super(WrapitLoader, self).__init__(mod_dict)
        self.args = args
        
    def _load_from_module(self, module, cmd_list):
        """load task from a module or dict with module members"""
        if inspect.ismodule(module):
            members = dict(inspect.getmembers(module))
        else:
            members = module
        task_list = self._load_tasks(members, cmd_list)
        doit_config = loader.load_doit_config(members)
        if doit_config == {}:
            doit_config = config_from_args(self.args)
        return task_list, doit_config

    def _load_tasks(self, members, cmd_list=()):
        """Get task generators and generate tasks

        @param members: (dict) containing the task generators, it might
                            contain other stuff
        @param command_names: (list - str) blacklist for task names
        @return task_list (list) of Tasks in the order they were defined on the file
        """

        # get functions defined in the module and select the task generators
        # a task generator function name starts with the string TASK_STRING
        funcs = []
        prefix_len = len(TASK_STRING)
        # get all functions defined in the module
        for name, ref in six.iteritems(members):

            # function is a task creator because of its name
            if inspect.isfunction(ref) and name.startswith(TASK_STRING):
                # remove TASK_STRING prefix from name
                task_name = name[prefix_len:]

            # object is a task creator because it contains the special method
            elif hasattr(ref, 'create_doit_tasks'):
                ref = ref.create_doit_tasks
                # If create_doit_tasks is a method, it should be called only
                # if it is bounded to an object.
                # This avoids calling it for the class definition.
                argspec = inspect.getargspec(ref)
                if len(argspec.args) != (1 if is_bound_method(ref) else 0):
                    continue
                task_name = name

            # ignore functions that are not a task creator
            elif True: # coverage can't get "else: continue"
                continue

            # tasks cant have name of commands
            if task_name in cmd_list:
                msg = ("Task can't be called '%s' because this is a command name."+
                       " Please choose another name.")
                raise InvalidDodoFile(msg % task_name)
            # get line number where function is defined
            line = inspect.getsourcelines(ref)[1]
            # add to list task generator functions
            funcs.append((task_name, ref, line))

        # sort by the order functions were defined (line number)
        # TODO: this ordering doesnt make sense when generators come
        # from different modules
        funcs.sort(key=lambda obj:obj[2])

        task_list = []
        for name, ref, line in funcs:
            try:
                task_return = ref(self.args)
            except TypeError:
                task_return = ref()
            generated_tasks = generate_tasks(name, task_return, ref.__doc__)
            for t in generated_tasks:
                t.options = vars(self.args)
            task_list.extend(generated_tasks)
        return task_list
