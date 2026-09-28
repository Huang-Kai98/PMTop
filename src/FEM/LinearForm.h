// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

//
// Created by huangkai on 25-1-21.
//

#ifndef LINEARFORM_H
#define LINEARFORM_H

#include "AnsysLoad.h"

class LinearForm {
public:
    LinearForm();

private:
    AnsysLoad fes;
};


#endif //LINEARFORM_H
