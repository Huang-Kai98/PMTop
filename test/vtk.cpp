//
// Created by huangkai on 25-1-21.
//

#include <AnsysLoad.h>
#include "Optimizer.h"

int main() {

    std::string datafile = "Mix_Shell_Plane2";
    AnsysLoad anasys(datafile);

    anasys.saveToVTK("Mix_Shell_Plane2.vtk");
    return 0;
}