#ifndef DERIVEDCONTAINER_H
#define DERIVEDCONTAINER_H

//template<typename U> class Container {
//public:
//  class iterator {
//    virtual U &operator *() = 0;
//  };
//};


//template<typename T, typename U> class DerivedContainer : public Container<U> {
template<typename T, typename U> class DerivedContainer {
public:
  typedef typename T::iterator TIt;
  DerivedContainer(T &container) : _container(container) {}
  
  class iterator {
//  class iterator : public Container<U>::iterator {
  public:
      iterator(TIt it) : _it(it) {}
//      virtual U &operator *() { return (U)(*_it); }
      U operator *() { return (U)(*_it); }
      TIt operator ++() { return ++_it; }
      bool operator !=(const iterator &o) { return _it!=o._it; }
  private:
      TIt _it;
  };

  iterator begin() { return iterator(_container.begin()); }
  iterator end() { return iterator(_container.end()); }

private:
  T &_container;
};


template<typename T, typename U> class ConstDerivedContainer {
public:
  typedef typename T::const_iterator TIt;
  ConstDerivedContainer(T &container) : _container(container) {}
  
  class const_iterator {
  public:
      const_iterator(TIt it) : _it(it) {}
      U operator *() { return (U)(*_it); }
      TIt operator ++() { return ++_it; }
      bool operator !=(const const_iterator &o) { return _it!=o._it; }
  private:
      TIt _it;
  };

  const_iterator begin() { return const_iterator(_container.begin()); }
  const_iterator end() { return const_iterator(_container.end()); }
  typename T::size_type size() { return _container.size(); }

private:
  T &_container;
};

#endif /* DERIVEDCONTAINER_H */

