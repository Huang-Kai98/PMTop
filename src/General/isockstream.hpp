// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once

#include <signal.h>
#include <sstream>
using namespace std;

/** Data type for input socket stream class. The class is used as server
    to receive data from a client on specified port number. The user gets
    data from the stream as from any other input stream. */
class isockstream {
private:
  int portnum, portID, socketID, error;
  char *Buf;

  int establish();
  int read_data(int socketid, char *buf, int size);

public:
  /** The constructor takes as input the portnumber port on which
      it establishes a server. */
  isockstream(int port);

  bool good() { return (!error); }

  /// Start waiting for data and return it in an input stream.
  void receive(istringstream **in);

  /** Virtual destructor. If the data hasn't been sent it sends it. */
  ~isockstream();
};