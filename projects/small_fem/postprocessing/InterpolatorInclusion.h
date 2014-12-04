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
                                       fullMatrix<scalar>& values,
                                       std::vector<bool>& isValid){
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

  // Alloc isValid //
  isValid.resize(nPoint);

  // Coef iterators //
           std::vector<Dof>                      dof;
  typename std::map<Dof, scalar>::const_iterator end = coef.end();
  typename std::map<Dof, scalar>::const_iterator it;

  // Iterate on 'point' //
  for(size_t i = 0; i < nPoint; i++){
    // Search all the elements containg this point
    SPoint3   thisPoint(point(i, 0), point(i, 1), point(i, 2));
    std::vector<MElement*> element =
      model.getMeshElementsByCoord(thisPoint, dim, true);

    // If no element found, point is invalid
    if(element.empty())
      isValid[i] = false;

    else{
      // Search an element of this GroupOfElement
      size_t nElement = element.size();
      size_t idx = 0;

      for(size_t e = 0; e < nElement; e++)
        if(goe.isMember(*element[idx]))
          break;
        else
          idx++;

      // Is point valid ?
      isValid[i] = !(idx == nElement);

      // If it is valid, proceed to interpolation
      if(isValid[i]){
        // Get Dofs related to this Element
        fs.getKeys(*element[idx], dof);
        const size_t size = dof.size();

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
            throw Exception("Interpolator: "
                            "missing Dof %s on element %d of GroupOfElement %d",
                            dof[j].toString().c_str(), i, goe.getId());
        }

        // Get Node coordinate
        fullVector<double> xyz(3);
        xyz(0) = point(i, 0);
        xyz(1) = point(i, 1);
        xyz(2) = point(i, 2);

        // Interpolate (AT LAST !!)
        fullVector<scalar> tmp;
        interpolate(*element[idx], fs, fsCoef, xyz, tmp);

        values(i, 0) = tmp(0);

        if(!isScalar){
          values(i, 1) = tmp(1);
          values(i, 2) = tmp(2);
        }
      }
    }
  }
}

template<typename scalar>
void Interpolator<scalar>::
interpolate(const GroupOfElement& goe,
            const GroupOfElement& point,
            const FunctionSpace& fs,
            const std::map<Dof, scalar>& coef,
            std::map<const MVertex*, std::vector<scalar> >& data){

  // Get the Vertices of 'point' //
  std::set<const MVertex*, VertexComparator> vertex;
  point.getAllVertex(vertex);

  // Get those Vertices coordinates //
  const size_t nVertex = vertex.size();
  fullMatrix<double> coordinate(nVertex, 3);

        std::set<const MVertex*, VertexComparator>::iterator it;
  const std::set<const MVertex*, VertexComparator>::iterator end = vertex.end();

  it = vertex.begin();
  for(size_t i = 0; it != end; it++, i++){
    coordinate(i, 0) = (*it)->x();
    coordinate(i, 1) = (*it)->y();
    coordinate(i, 2) = (*it)->z();
  }

  // Interpolate //
  fullMatrix<scalar>  value;
  std::vector<bool> isValid;
  interpolate(goe, fs, coef, coordinate, value, isValid);

  // Get Data Map //
  const size_t nDim = value.size2();
  std::vector<scalar> tmp(nDim);

  it = vertex.begin();
  for(size_t i = 0; it != end; it++, i++){
    // Insert Vertex only if valid
    if(isValid[i]){
      // Pair
      std::pair<const MVertex*, std::vector<scalar> > pair;

      // Vertex
      pair.first = (*it);

      // Value
      pair.second.resize(nDim);
      for(size_t j = 0; j < nDim; j++)
        pair.second.at(j) = value(i, j);

      // Add to data
      data.insert(pair);
    }
  }
}

template<typename scalar>
void Interpolator<scalar>::dump(std::ofstream& stream,
                                const std::map<int,std::vector<double> >& data){
  // Iterator //
  std::map<int, std::vector<double> >::const_iterator  it = data.begin();
  std::map<int, std::vector<double> >::const_iterator end = data.end();

  // Dump //
  for(; it != end; it++){
    stream << it->first << " ";

    for(size_t i = 0; i < it->second.size(); i++)
      stream << it->second[i] << " ";

    stream << std::endl;
  }
}
