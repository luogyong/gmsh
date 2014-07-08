#ifndef _MPIOSTREAM_H_
#define _MPIOSTREAM_H_

#include <ostream>

/**
   @class MPIOStream
   @brief OStream for MPI

   This class handles OStreams for multiple MPI process:
   only a given rank can access to it.
*/

class MPIOStream{
 private:
  bool myStream;
  std::ostream* ostream;

 public:
   MPIOStream(int rank, std::ostream& ostream);
  ~MPIOStream(void);

  template <typename T>
    const MPIOStream& operator<<(const T& t) const;

  const MPIOStream&   operator<<(std::ostream& (*t)(std::ostream&)) const;
};

/**
   @fn MPIOStream::MPIOStream
   @param rank The MPI process rank that will be allowed to use this MPIStream
   @param ostream The OStream this MPI process will use

   Instanciates a new MPIOStream
   **

   @fn MPIOStream::~MPIOStream

   Deltes this MPIOStream
   **

   @fn MPIOStream::operator<<(const T& t) const
   @param t The object or value to be inserted in this MPIOStream
   @return Returns this MPIOStream

   Inserts the given argument into this MPIOStream
   **

   @fn MPIOStream::operator<<(std::ostream& (*t)(std::ostream&)) const
   @param t A ostream function to be inserted in this MPIOStream
   @return Returns this MPIOStream

   Inserts the given argument into this MPIOStream
 */

/////////////////////////////////
// Template & Inline Functions //
/////////////////////////////////
template <typename T>
inline const MPIOStream& MPIOStream::operator<<(const T& t) const{
  if(myStream)
    *ostream << t;

  return *this;
}

inline const MPIOStream& MPIOStream::
operator<<(std::ostream& (*t)(std::ostream&)) const{
  if(myStream)
    *ostream << t;

  return *this;
}

#endif
