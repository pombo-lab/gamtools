import sys

from .loader import WrapitLoader
from doit.doit_cmd import DoitMain

def run(task_creators, args, task_selectors=[]):
    """run doit using task_creators

    @param task_creators: module or dict containing task creators
    """
    sys.exit(DoitMain(WrapitLoader(args, task_creators)).run(task_selectors))
