// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "osockstream.hpp"
#include <iostream>
#include <netdb.h>
#include <netinet/in.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

osockstream::osockstream(int port, const char *hostname) : ostringstream() {
  int n;

  sent = 0;
  portnum = port;
  n = strlen(hostname);
  if (n >= 128) {
    strncpy(machine, hostname, 127);
    machine[127] = '\0';
  } else
    strcpy(machine, hostname);
}

int osockstream::call_socket() {
  struct sockaddr_in sa;
  struct hostent *hp;
  int s;

  if ((hp = gethostbyname(machine)) == NULL)
    return (-1);

  memset(&sa, 0, sizeof(sa));
  memcpy((char *)&sa.sin_addr, hp->h_addr, hp->h_length); /* set address */
  sa.sin_family = hp->h_addrtype;
  sa.sin_port = htons((u_short)portnum);

  if ((s = socket(hp->h_addrtype, SOCK_STREAM, 0)) < 0) /* get socket */
    return (-2);
  if (connect(s, (const sockaddr *)&sa, sizeof sa) < 0) { /* connect */
    close(s);
    return (-3);
  }
  return (s);
}

int osockstream::send_through_socket(int s, const char *Buf, int size) {
  int bcount = 0; /* counts bytes read */
  int br = 0;     /* bytes read this pass */
  const char *b = Buf;

  char length[32];
  sprintf(length, "%d", size);

  ::write(s, length, 32);

  while (bcount < size)
    if ((br = ::write(s, b, size - bcount)) > 0) {
      bcount += br;
      b += br;
    } else if (br < 0)
      return (-1);
  return 1;
}

int osockstream::send() {
  int sock, size, ierr = 0;
  const char *Buf;

  if (sent)
    return (-1); // Data already sent
  sent = 1;
  switch (sock = call_socket()) {
  case -1:
    cerr << "Unknown host: " << machine << endl;
    ierr = -2;
    break;
  case -2:
    cerr << "Unable to establish socket on the local machine" << endl;
    ierr = -3;
    break;
  case -3:
    cerr << "Unable to connect to port " << portnum << " on " << machine
         << endl;
    ierr = -4;
    break;
  default: {
    string my_str(str());
    Buf = my_str.c_str();
    size = my_str.size();
    switch (send_through_socket(sock, Buf, size)) {
    case -1:
      cerr << "Error sending data to port " << portnum << " on " << machine
           << endl;
      ierr = -5;
      break;
    default:
      break;
    }
    close(sock);
  } break;
  }
  return ierr;
}

osockstream::~osockstream() {
  if (!sent)
    send();
}
