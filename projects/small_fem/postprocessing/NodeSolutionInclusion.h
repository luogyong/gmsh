/////////////////////////////////////////////////
// Templates Implementations for NodeSolution: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

template<typename scalar>
NodeSolution<scalar>::NodeSolution(void){
  pView = new PViewDataGModel(PViewDataGModel::NodeData);
}

template<typename scalar>
NodeSolution<scalar>::~NodeSolution(void){
  pView->destroyData();
  delete pView;
}

template<typename scalar>
void NodeSolution<scalar>::write(std::string fileName) const{
  pView->setName(fileName);
  pView->writeMSH(fileName + ".msh");
}
