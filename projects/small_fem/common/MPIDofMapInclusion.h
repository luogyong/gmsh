//////////////////////////////////////////////
// Templates Implementations for MPIDofMap: //
// Inclusion compilation model              //
//                                          //
// Damn you gcc: we want 'export' !         //
//////////////////////////////////////////////

template<typename T>
MPIDofMap<T>::MPIDofMap(void){
}

template<typename T>
MPIDofMap<T>::~MPIDofMap(void){
}

template<typename T>
void MPIDofMap<T>::
getGlobalMap(const std::map<Dof, T>& local, std::map<Dof, T>& global){
  // Serialize Dofs //
  int* entity;
  int* type;
  int  size = serialize(local, &entity, &type);

  // Gather sizes from other Nodes //
  int  sizeSum;
  int* sizeAll = gatherSize(size, &sizeSum);

  // Strides //
  int* stride = cpteStride(sizeAll);

  // Exchange //
  int* allEntity = exchange(entity, size, sizeSum, sizeAll, stride);
  int* allType   = exchange(type,   size, sizeSum, sizeAll, stride);

  // Create Map //
  unserialize(global, allEntity, allType, sizeSum);

  // Clear //
  delete[] entity;
  delete[] type;
  delete[] sizeAll;
  delete[] stride;
  delete[] allEntity;
  delete[] allType;
}

template<typename T>
void MPIDofMap<T>::
getDofOwners(const std::map<Dof, T>& local, std::multimap<Dof, int>& owners){
  // Serialize Dofs //
  int* entity;
  int* type;
  int  size = serialize(local, &entity, &type);

  // Local owner vector //
  int  myProc;
  int* owner = new int[size];

  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  for(int i = 0; i < size; i++)
    owner[i] = myProc;

  // Gather sizes from other Nodes //
  int  sizeSum;
  int* sizeAll = gatherSize(size, &sizeSum);

  // Strides //
  int* stride = cpteStride(sizeAll);

  // Exchange //
  int* allEntity = exchange(entity, size, sizeSum, sizeAll, stride);
  int* allType   = exchange(type,   size, sizeSum, sizeAll, stride);
  int* allOwners = exchange(owner,  size, sizeSum, sizeAll, stride);

  // Create Map //
  unserialize(owners, allEntity, allType, allOwners, sizeSum);

  // Clear //
  delete[] entity;
  delete[] type;
  delete[] owner;
  delete[] sizeAll;
  delete[] stride;
  delete[] allEntity;
  delete[] allType;
  //delete[] allOwners;
}

template<typename T>
int MPIDofMap<T>::
serialize(const std::map<Dof, T>& in, int** entity, int** type){
  // Get Size //
  int size = in.size();

  // Allocate //
  *entity = new int[size];
  *type   = new int[size];

  // Serialize //
  typename std::map<Dof, T>::const_iterator it;
  typename std::map<Dof, T>::const_iterator end;

  // Entity
  it  = in.begin();
  end = in.end();
  for(int i = 0; it != end; it++, i++)
    (*entity)[i] = it->first.getEntity();

  // Type
  it  = in.begin();
  end = in.end();
  for(int i = 0; it != end; it++, i++)
    (*type)[i] = it->first.getType();

  // Done //
  return size;
}

template<typename T>
void MPIDofMap<T>::
unserialize(std::map<Dof, T>& map, int* entity, int* type, int size){
  // Clear //
  map.clear();

  // Build //
  for(int i = 0; i < size; i++)
    map.insert(std::pair<Dof, T>(Dof(entity[i], type[i]), 0));
}

template<typename T>
void MPIDofMap<T>::unserialize(std::multimap<Dof, int>& map,
                               int* entity, int* type, int* owner, int size){
  // Clear //
  map.clear();

  // Build //
  for(int i = 0; i < size; i++)
    map.insert(std::pair<Dof, int>(Dof(entity[i], type[i]), owner[i]));
}

template<typename T>
int* MPIDofMap<T>::gatherSize(int mySize, int* sum){
  // Number of MPI Process //
  int nProc;
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  // Alloc //
  int* size = new int[nProc];

  // Gather //
  MPI_Allgather(&mySize, 1, MPI_INT, size, 1, MPI_INT, MPI_COMM_WORLD);

  // Sum //
  *sum = 0;
  for(int i = 0; i < nProc; i++)
    *sum += size[i];

  // Done //
  return size;
}

template<typename T>
int* MPIDofMap<T>::cpteStride(int* size){
  // Number of MPI Process //
  int nProc;
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  // Alloc //
  int* stride = new int[nProc];

  // Populate //
  stride[0] = 0;
  for(int i = 1; i < nProc; i++)
    stride[i] = stride[i - 1] + size[i - 1];

  // Done //
  return stride;
}

template<typename T>
int* MPIDofMap<T>::
exchange(int* myData, int mySize, int sizeSum, int* size, int* stride){
  // Alloc //
  int* allData = new int[sizeSum];

  // All Gather //
  MPI_Allgatherv(myData , mySize,         MPI_INT,
                 allData,   size, stride, MPI_INT,
                 MPI_COMM_WORLD);
  // Done //
  return allData;
}
