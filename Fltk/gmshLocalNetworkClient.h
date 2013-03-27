// Gmsh - Copyright (C) 1997-2013 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.

#ifndef _GMSH_LOCAL_NETWORK_CLIENT_H_
#define _GMSH_LOCAL_NETWORK_CLIENT_H_

#include <vector>
#include <algorithm>
#include "onelab.h"

class gmshLocalNetworkClient : public onelab::localNetworkClient{
 private:
  // a gmsh local network client can launch subclients (this is typical for a
  // metamodel that calls several underlying models); _clients keeps track of
  // the master (this) and the subclients.
  std::vector<gmshLocalNetworkClient*> _clients;
 public:
  gmshLocalNetworkClient(const std::string &name, const std::string &executable,
                         const std::string &remoteLogin="")
    : onelab::localNetworkClient(name, executable, remoteLogin)
  {
    addClient(this);
  }
  void addClient(gmshLocalNetworkClient *client)
  {
    _clients.push_back(client);
  }
  void removeClient(gmshLocalNetworkClient *client)
  {
    std::vector<gmshLocalNetworkClient*>::iterator it;
    it = std::find(_clients.begin(), _clients.end(), client);
    if(it != _clients.end()) _clients.erase(it);
  }
  int getNumClients(){ return _clients.size(); }
  gmshLocalNetworkClient *getClient(int i)
  {
    if(i >= 0 && i < getNumClients()) return _clients[i];
    return 0;
  }
  bool receiveMessage(int &type);
  bool run();
  bool kill();
};

#endif
