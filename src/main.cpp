#include <iostream>
#include <vector>
#include "reflow.hpp"

int main(int argc, char** argv)
{
    reflow S(20,3,std::vector<double>{1,10,100000});
    S.solve();

    return 0;
}
