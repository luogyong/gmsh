/////////////////////////////////////////////////
// Templates Implementations for Interpolator: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

template<typename scalar>
Interpolator<scalar>::Interpolator(void){
}

template<typename scalar>
Interpolator<scalar>::~Interpolator(void){
}

template<typename scalar>
void Interpolator<scalar>::interpolate(const GroupOfElement& goe,
                                       const FunctionSpace& fs,
                                       const std::map<Dof, scalar>& coef,
                                       const fullMatrix<double>& point,
                                       fullMatrix<scalar>& values){
  // Get GModel //
  GModel&    model = goe.getMesh().getModel();
  const size_t dim = model.getDim();

  // Scalar or Vector ?
  const bool isScalar = fs.isScalar();

  // Alloc values //
  const size_t nPoint = point.size1();

  if(isScalar)
    values.resize(nPoint, 1);
  else
    values.resize(nPoint, 3);

  // Coef iterators //
  const typename std::map<Dof, scalar>::const_iterator end = coef.end();
  typename std::map<Dof, scalar>::const_iterator it;

  // Iterate on 'point' //
  for(size_t i = 0; i < nPoint; i++){
    // Search element containg this point
    SPoint3   thisPoint(point(i, 0), point(i, 1), point(i, 2));
    MElement* element = model.getMeshElementByCoord(thisPoint, dim, true);

    // WARNING: if no element found, set 'values' to zero
    // Note: 'zero' must have a meaning for the templatet scalar type
    if(!element){
      values(i, 0) = 0;

      if(!isScalar){
        values(i, 1) = 0;
        values(i, 2) = 0;
      }
    }

    else{
      // Get Dofs related to this Element
      const std::vector<Dof>& dof  = fs.getKeys(*element);
      const size_t            size = dof.size();

      // Get Coef In FS Basis
      std::vector<scalar> fsCoef(size);
      for(size_t j = 0; j < size; j++){
        // Get Value of Dof 'j'
        it = coef.find(dof[j]);

        // If found in map
        if(it != end)
          fsCoef[j] = it->second;

        // Else
        else
          fsCoef[j] = 0;
      }

      // Get Node coordinate
      fullVector<double> xyz(3);
      xyz(0) = point(i, 0);
      xyz(1) = point(i, 1);
      xyz(2) = point(i, 2);

      // Interpolate (AT LAST !!)
      fullVector<scalar> tmp;
      interpolate(*element, fs, fsCoef, xyz, tmp);

      values(i, 0) = tmp(0);

      if(!isScalar){
        values(i, 1) = tmp(1);
        values(i, 2) = tmp(2);
      }
    }
  }
}
