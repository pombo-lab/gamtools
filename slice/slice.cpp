#include <iostream>
#include <string>
using namespace std;

string my_function(string prefix)
{
    string output;
    output = prefix + "balls";
    return output;
}

int main()
{
    string a;
    a = my_function("Mr. ");
    cout << a;
    cout << "\n";
    return 0;
}
