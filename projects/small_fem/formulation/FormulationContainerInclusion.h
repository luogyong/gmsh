/////////////////////////////////////////////////////////
// Templates Implementations for FormulationContainer: //
// Inclusion compilation model                         //
//                                                     //
// Damn you gcc: we want 'export' !                    //
/////////////////////////////////////////////////////////

template<typename scalar>
FormulationContainer<scalar>::FormulationContainer(void){
}

template<typename scalar>
FormulationContainer<scalar>::~FormulationContainer(void){
}

template<typename scalar>
void FormulationContainer<scalar>::
addFormulation(FormulationBlock<scalar>& formulation){
  fList.push_back(&formulation);
}

template<typename scalar>
void FormulationContainer<scalar>::
addFormulation(FormulationCoupled<scalar>& formulation){
  const std::list<FormulationBlock<scalar>*>
    tmp = formulation.getFormulationBlocks();

  typename std::list<FormulationBlock<scalar>*>::const_iterator
    it  = tmp.begin();
  typename std::list<FormulationBlock<scalar>*>::const_iterator
    end = tmp.end();

  for(; it != end; it++)
    fList.push_back(*it);
}

template<typename scalar>
const std::list<FormulationBlock<scalar>*>&
FormulationContainer<scalar>::getFormulationBlocks(void) const{
  return fList;
}

template<typename scalar>
void FormulationContainer<scalar>::update(void){
  typename std::list<FormulationBlock<scalar>*>::iterator
    it  = fList.begin();
  typename std::list<FormulationBlock<scalar>*>::iterator
    end = fList.end();

  for(; it != end; it++)
      (*it)->update();
}
