#ifndef MHS_VECTOR_H
#define MHS_VECTOR_H

#ifndef M_PI
  #define M_PI		3.14159265358979323846
#endif  

#include <stdio.h>
#include <string.h>
#include "rand_help.h"
#include <sstream>



template <typename T>
T get_dist_from_line(const Vector<T,3> & lx1,const Vector<T,3> & lx2,const Vector<T,3> & x3,T max_dist=-1)
{
  Vector<T,3> n;
  n=lx2-lx1;
  T l_n=std::sqrt(n.norm2())+std::numeric_limits< T >::epsilon();
  n/=l_n;
  
  Vector<T,3> v_lx1_x3;
  v_lx1_x3=x3-lx1;
  
  T l_v1=v_lx1_x3.dot(n);
  
  T d;
  
  if (l_v1>0)
  {
   if (l_v1<l_n)
   {
      d=std::sqrt((v_lx1_x3-n*l_v1).norm2());
      
   }else
   {
      d=std::sqrt((x3-lx2).norm2());
   }
  }else
  {
      d=std::sqrt((x3-lx1).norm2());
  }
  
  if (max_dist>0)
  {
    return  std::min(d,max_dist);
  }
  
  return d;
}



//NOTE distance weighted with thickness (0 at centerline, 1 at surface)
template <typename T>
T get_dist_from_surface(const Vector<T,3> & lx1,T s1,const Vector<T,3> & lx2,T s2,const Vector<T,3> & x3,T max_dist=-1)
{
  Vector<T,3> n;
  n=lx2-lx1;
  T l_n=std::sqrt(n.norm2())+std::numeric_limits< T >::epsilon();
  n/=l_n;
  
  Vector<T,3> v_lx1_x3;
  v_lx1_x3=x3-lx1;
  
  T l_v1=v_lx1_x3.dot(n);
  
  T d;
  
  if (l_v1>0)
  {
   if (l_v1<l_n)
   {
      d=std::sqrt((v_lx1_x3-n*l_v1).norm2());
      T w=(l_n-l_v1)/l_n;
      d=d/(w*s1+(1-w)*s2);
      
   }else
   {
      d=std::sqrt((x3-lx2).norm2())/s2;
   }
  }else
  {
      d=std::sqrt((x3-lx1).norm2())/s1;
  }
  
  if (max_dist>0)
  {
    return  std::min(d,max_dist);
  }
  
  return d;
}


template<typename T,int Dim>
class Vector
{
public:
    T v[Dim];
    Vector() {};
    
    
    template<int Dim2>
    Vector<T,Dim2> get_subvec()
    {
      Vector<T,Dim2> vec;
      for (int i=0; i<std::min(Dim,Dim2); i++)
	vec[i]=v[i];
      return vec;
	
    }

    void print() const
    {
        std::stringstream s;
        s<<"(";

        for (int i=0; i<Dim; i++)
            if (i<Dim-1)
                s<<v[i]<<",";
            else
                s<<v[i]<<")";
        printf("%s\n",s.str().c_str());
    }

    const T &operator[](int i) const {
        return v[i];
    };
    T &operator[](int i)  {
        return v[i];
    };

    //const T &operator[](int i)  { return v[i]; };

    T norm2() const
    {
        T n=0;
        #pragma omp simd reduction(+:n)
        for (int i=0; i<Dim; i++)
            n+=v[i]*v[i];
        return n;
    }

    T norm1() const
    {
        T n=0;
        #pragma omp simd reduction(+:n)
        for (int i=0; i<Dim; i++)
            n+=std::abs(v[i]);
        return n;
    }

    void normalize()
    {
        (*this )/=(std::sqrt(this->norm2())+std::numeric_limits<T>::epsilon());
    }

    const Vector & rand(T min,T max)
    {
        for (int i=0; i<Dim; i++)
            v[i]=myrand(max-min)+min;
        return *this;
    }

    const Vector & rand(T* min,T* max)
    {
        for (int i=0; i<Dim; i++)
            v[i]=myrand(max[i]-min[i])+min[i];
        return *this;
    }

    const Vector & rand_normal(T sigma)
    {
        for (int i=0; i<Dim; i++)
            v[i]=mynormalrand(sigma);
        return *this;
    }
    
    const Vector & rand_normal(Vector sigma)
    {
        for (int i=0; i<Dim; i++)
            v[i]=mynormalrand(sigma[i]);
        return *this;
    }

    Vector(const T rhs[]) {
        memcpy (v,rhs,Dim*sizeof(T));
    };



    Vector(const T z,const T y=0,const T x=0, const T w=0)
    {
        v[0]=z;
        if (Dim>1)
            v[1]=y;
        if (Dim>2)
            v[2]=x;
	if (Dim>3)
            v[3]=w;
    };
    
    
    
    
    Vector & set(const T z,const T y=0,const T x=0)
    {
	v[0]=z;
        if (Dim>1)
            v[1]=y;
        if (Dim>2)
            v[2]=x;
        return *this;
    }    

    
    bool all_greater(const Vector &rhs) const
    {
        bool result=true;
        for (int i=0; i<Dim; i++)
	{
	    result&=v[i]>rhs[i];
	}
        return result;
    }
    
    bool all_smaller(const Vector &rhs) const
    {
        bool result=true;
        for (int i=0; i<Dim; i++)
	{
	    result&=v[i]<rhs[i];
	}
        return result;
    }
    
    bool all_unequal(const Vector &rhs) const
    {
        bool result=true;
        for (int i=0; i<Dim; i++)
	{
	    result&=(v[i]!=rhs[i]);
	}
        return result;
    }
    
    bool any_greater(const Vector &rhs) const
    {
        bool result=false;
        for (int i=0; i<Dim; i++)
	{
	    result|=v[i]>rhs[i];
	}
        return result;
    }
    
    bool any_smaller(const Vector &rhs) const
    {
        bool result=false;
        for (int i=0; i<Dim; i++)
	{
	    result|=v[i]<rhs[i];
	}
        return result;
    }
    
    bool any_unequal(const Vector &rhs) const
    {
        bool result=false;
        for (int i=0; i<Dim; i++)
	{
	    result|=(v[i]!=rhs[i]);
	}
        return result;
    }
    
    

    T dot(const Vector &rhs) const
    {
        T result=0;
        #pragma omp simd reduction(+:result)
        for (int i=0; i<Dim; i++)
            result+=v[i]*rhs.v[i];
        return result;
    }

    template<typename S>
    T dot(const S rhs[]) const
    {
        T result=0;
        #pragma omp simd reduction(+:result)
        for (int i=0; i<Dim; i++)
            result+=v[i]*rhs[i];
        return result;
    }

    const Vector cross(const Vector &rhs) const
    {
        sta_assert_error(Dim==3);
        Vector c;
        c.v[0]=v[1]*rhs.v[2]-v[2]*rhs.v[1];
        c.v[1]=v[2]*rhs.v[0]-v[0]*rhs.v[2];
        c.v[2]=v[0]*rhs.v[1]-v[1]*rhs.v[0];
        return c;
    }

    Vector & operator=(const Vector &rhs)
    {
        if (this != &rhs)
            memcpy (v,rhs.v,Dim*sizeof(T));
        return *this;
    }
    
    template<typename R>
    Vector & operator=(const Vector<R,Dim> &rhs)
    {
        for (int i=0;i<3;i++)
	  v[i]=rhs.v[i];
        return *this;
    }    

    template<typename S>
    Vector & operator=(const std::vector<S> &rhs)
    {
        for (int i=0; i<Dim; i++)
        {
            if (i<rhs.size())
                v[i]=rhs[i];
            else
                v[i]=0;
        }
        return *this;
    }

    Vector & operator=(const T rhs[])
    {
        if (v != rhs)
            //    if (this != &rhs)
            memcpy (v,rhs,Dim*sizeof(T));
        return *this;
    }


//     template<typename S>
//     Vector & operator=(S f)
//     {
//         #pragma omp simd
//         for (int i=0; i<Dim; i++)
//             v[i]=f;
//         return *this;
//     }
    
    Vector & operator=(T f)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]=f;
        return *this;
    }

    template<typename D>
    Vector & operator+=(const Vector<D,Dim> &rhs)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]+=rhs.v[i];
        return *this;
    }

    Vector & operator+=(const T rhs[])
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]+=rhs[i];
        return *this;
    }

    Vector & operator+=(T f)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]+=f;
        return *this;
    }

    Vector & operator-=(const Vector &rhs)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]-=rhs.v[i];
        return *this;
    }

    Vector & operator-=(T f)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]-=f;
        return *this;
    }

    Vector & operator*=(const Vector &rhs)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]*=rhs.v[i];
        return *this;
    }

    Vector & operator*=(T f)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]*=f;
        return *this;
    }

    template<typename D>
    Vector & operator/=(const Vector<D,Dim> &rhs)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]/=rhs.v[i];
        return *this;
    }

    Vector & operator/=(T f)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]/=f;
        return *this;
    }

    template<typename D>
    const Vector operator+(const Vector<D,Dim> &rhs) const
    {
        return Vector(*this) += rhs;
    }
    
    

    const Vector operator+(const T rhs[]) const
    {
        return Vector(*this) += rhs;
    }

    const Vector operator+(T f) const
    {
        return Vector(*this) += f;
    }

    const Vector operator-(const Vector &rhs) const
    {
        return Vector(*this) -= rhs;
    }

    const Vector operator-(T f) const
    {
        return Vector(*this) -= f;
    }

     template<typename D>
    const Vector operator/(const Vector<D,Dim> &rhs) const
    {
        return Vector(*this) /= rhs;
    }

    const Vector operator/(T f) const
    {
        return Vector(*this) /= f;
    }

    const Vector operator*(const Vector &rhs) const
    {
        return Vector(*this) *= rhs;
    }

    const Vector operator*(T f) const
    {
        return Vector(*this) *= f;
    }
    
    Vector & max(const Vector & v2)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]=std::max(v[i],v2.v[i]);
        return *this;
    }
    
    Vector & min(const Vector & v2)
    {
        #pragma omp simd
        for (int i=0; i<Dim; i++)
            v[i]=std::min(v[i],v2.v[i]);
        return *this;
    }

    bool operator==(const Vector &other) const {
        for (int i=0; i<Dim; i++)
            if (v[i]!=other.v[i])
                return false;
        return true;
    }

    bool operator!=(const Vector &other) const {
        return !(*this==other);
    }
};









#endif
