#include <iostream>
#include "utils.h"
using namespace std;

int main(const int argc, const char** argv) {
    cout.precision(6);
    cout << "hmp" << endl;
    clock_t at, ad;

    cout << "Parse Parameters" << endl;
    if (argc == 1)
    {
        helpmsg(argv[0]);
        return -1;
    }

    const int k = atoi(read(argc, argv, "-k", ""));
    int dim = atoi(read(argc, argv, "-d", ""));
    const char* methodName = read(argc, argv, "-m", "");
    const char* datafile = read(argc, argv, "-f", "");  // option file name
    int h = atoi(read(argc, argv, "-h", ""));



    return 0;
}
