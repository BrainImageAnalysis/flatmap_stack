#ifndef RAND_HELP_H
#define RAND_HELP_H
#include <complex>
#include <cmath>
#include "mhs_vector.h"
#include <random>


template<typename T,int Dim> class Vector;

int maxpiter=0;
std::default_random_engine generator;
std::normal_distribution<double> normal_distribution(0.0,1);


template<typename T,typename TData,int Dim>
bool interp3D( T & result,
              TData * image,
              const Vector<T,Dim>& position,
              const std::size_t shape[])
{
    sta_assert_debug0(image!=NULL);
    std::size_t  indexes[8];
    T weights[8];

    if (!sub2indInter(
                shape,
                position,
                indexes,
                weights))
    {
        return false;
    }
    result=0;
    for (int a=0; a<8; a++)
    {
        result+=weights[a]*image[indexes[a]];
    }
    return true;
}


template<typename D,typename T,int Dim>
std::size_t sub2ind(D shape[],Vector<T,Dim> & pos)
{
    return (std::floor(pos[0])*shape[1]
            +std::floor(pos[1]))*shape[2]
           +std::floor(pos[2]);
}

template<typename D,typename T,int Dim>
bool sub2indInter(const D shape[],
                  const Vector<T,Dim> & pos,
                  std::size_t  indexes[],
                  T weights[])
{
    int Z=std::floor(pos[0]);
    int Y=std::floor(pos[1]);
    int X=std::floor(pos[2]);

//   sta_assert_error(Z+1<shape[0]);
//   sta_assert_error(Y+1<shape[1]);
//   sta_assert_error(X+1<shape[2]);

    if ((Z+1>=shape[0])||(Y+1>=shape[1])||(X+1>=shape[2]))
        return false;

    T wz=(pos[0]-Z);
    T wy=(pos[1]-Y);
    T wx=(pos[2]-X);

    indexes[0]=(Z*shape[1]+Y)*shape[2]+X;
    indexes[1]=(Z*shape[1]+Y)*shape[2]+(X+1);
    indexes[2]=(Z*shape[1]+(Y+1))*shape[2]+X;
    indexes[3]=(Z*shape[1]+(Y+1))*shape[2]+(X+1);
    indexes[4]=((Z+1)*shape[1]+Y)*shape[2]+X;
    indexes[5]=((Z+1)*shape[1]+Y)*shape[2]+(X+1);
    indexes[6]=((Z+1)*shape[1]+(Y+1))*shape[2]+X;
    indexes[7]=((Z+1)*shape[1]+(Y+1))*shape[2]+(X+1);

    weights[0]=(1-wz)*(1-wy)*(1-wx);
    weights[1]=(1-wz)*(1-wy)*(wx);
    weights[2]=(1-wz)*(wy)*(1-wx);
    weights[3]=(1-wz)*(wy)*(wx);
    weights[4]=(wz)*(1-wy)*(1-wx);
    weights[5]=(wz)*(1-wy)*(wx);
    weights[6]=(wz)*(wy)*(1-wx);
    weights[7]=(wz)*(wy)*(wx);

    return true;
}





inline float myrand(double scale)
{
    return (scale*((double)std::rand()/(double)RAND_MAX));
}



inline float myrand(float min, float max)
{
    //return std::max(min,std::min(max,(min+myrand(max-min))));
  return (min+myrand(max-min));
}




// double rand_normal(double mean, double stddev) {
//     static double n2 = 0.0;
//     static int n2_cached = 0;
//     if (!n2_cached) {
//         double x, y, r;
//         do {
//             x = 2.0*rand()/RAND_MAX - 1;
//             y = 2.0*rand()/RAND_MAX - 1;
//             r = x*x + y*y;
//         } while (r == 0.0 || r > 1.0);
//         {
//             double d = sqrt(-2.0*log(r)/r);
//             double n1 = x*d;
//             n2 = y*d;
//             double result = n1*stddev + mean;
//             n2_cached = 1;
//             return result;
//         }
//     } else {
//         n2_cached = 0;
//         return n2*stddev + mean;
//     }
// }

template<typename T>
T mynormalrand(T sigma)
{
    //return rand_normal(0,sigma);
    return sigma*normal_distribution(generator);
}

template<typename T>
T mynormalrand(T mean, T sigma)
{
    //return rand_normal(mean,sigma);
    return sigma*normal_distribution(generator)+mean;
}



/*!
 *  randomly pic point position
 */
template<typename T,typename P,typename TData>
T rand_pic_position(TData * img,const std::size_t shape[],P pos[],std::size_t  &center,int stride=1,int offset=0,bool subpixel=false,T maxvalue=1)
{
    T choice=myrand(maxvalue);
    std::size_t left_pos=0;
    std::size_t numvoxel=shape[0]*shape[1]*shape[2];
    std::size_t right_pos=(numvoxel);
    center=right_pos;
    int it=0;

    while (left_pos+1<right_pos)
    {
        center=(int)std::ceil(((T)right_pos+(T)left_pos)/(T)2);
        if (img[center*stride+offset-1]>=choice)
            right_pos=center;
        else
            left_pos=center;
        it++;
    }
    center=right_pos-1;


    if (maxpiter<it)
        maxpiter=it;

    std::size_t indx=center;

    std::size_t posI[3];
    posI[0]=indx/(shape[1]*shape[2]);
    indx%=(shape[1]*shape[2]);
    posI[1]=indx/(shape[2]);
    posI[2]=indx%(shape[2]);
    if (subpixel)
    {
//         pos[0]=posI[0]+myrand(1-std::numeric_limits<T>::min());
//         pos[1]=posI[1]+myrand(1-std::numeric_limits<T>::min());
//         pos[2]=posI[2]+myrand(1-std::numeric_limits<T>::min());
	pos[0]=posI[0]+myrand(1-std::numeric_limits<T>::epsilon());
        pos[1]=posI[1]+myrand(1-std::numeric_limits<T>::epsilon());
        pos[2]=posI[2]+myrand(1-std::numeric_limits<T>::epsilon());
    } else
    {
        pos[0]=posI[0]+0.5;
        pos[1]=posI[1]+0.5;
        pos[2]=posI[2]+0.5;
    }
    return choice;
}



/*
 *
 *  randomly pic array pos
 *
 * */
//TODO maxvalue=img[num_elements-1];
template<typename T,typename TData>
std::size_t rand_pic_array(TData * img,std::size_t num_elements,T maxvalue=1)
{
    T choice=myrand(maxvalue);
    std::size_t left_pos=0;
    std::size_t right_pos=(num_elements);
    std::size_t center=right_pos;
    while (left_pos+1<right_pos)
    {
        center=(int)std::ceil(((T)right_pos+(T)left_pos)/(T)2);
        if (img[center-1]>=choice)
            right_pos=center;
        else
            left_pos=center;
    }
    return right_pos-1;
}

template<typename T,int Dim>
int random_pic_direction(Vector<T,Dim> & old_dir,
                         Vector<T,Dim> & new_dir,
                         T   sigma=1)
{

    Vector<T,Dim> rand_dir;
    rand_dir.rand_normal(sigma);
    new_dir=old_dir+rand_dir;
    new_dir.normalize();
    return -1;
}

template<typename T,int Dim>
 Vector<T,Dim> random_pic_direction(const Vector<T,Dim> & old_dir,
                         T   sigma=1)
{

    Vector<T,Dim> rand_dir;
    rand_dir.rand_normal(sigma);
    rand_dir+=old_dir;
    rand_dir.normalize();
    return rand_dir;
}


template<typename T,typename TData,int Dim>
int random_pic_direction(T n[],
                         int current=-1,
                         T   jump=0.25,
                         TData * sphere_pts=NULL,
                         TData * sphere_NN=NULL,
                         int num_sphere_pts=0)
{
    const int numdir=360;
    static T directions[2][numdir];
    static bool doinit=true;
    int index=-1;

    if (Dim==2)
    {

        if (doinit)
        {
            doinit=false;
            std::complex<T> w=exp(std::complex<T>(0,2*M_PI)/(T)(numdir));
            for (int i=0; i<numdir; i++)
            {
                std::complex<T> w1=std::pow(w,i);
                //printf("%f %f \n",w1.real(),w1.imag());
                directions[0][i]=w1.real();
                directions[1][i]=w1.imag();
            }
        }

        if (current==-1)
            index=std::rand()%numdir;
        else
        {
            int next=std::floor(jump*numdir);
            //printf("%f %d\n",jump,next);
            if (std::rand()%100>=50)
                index=(numdir+current+std::rand()%next)%numdir;
            else
                index=(numdir+current-std::rand()%next)%numdir;
        }

        sta_assert_error(index>=0);
        sta_assert_error(index<numdir);
        n[0]=directions[0][index];
        n[1]=directions[1][index];
    }

    if (Dim==3)
    {



        if (doinit)
        {
            doinit=false;
            sta_assert_error(sphere_pts!=NULL);
            sta_assert_error(sphere_NN!=NULL);
            sta_assert_error(num_sphere_pts>0);
        }

        if (current==-1)
            index=std::rand()%num_sphere_pts;
        else
        {
            //int next=std::max(std::floor(jump*num_sphere_pts),T(25));
            int next=std::max(std::floor(jump*num_sphere_pts),T(25));
            sta_assert_error(next<num_sphere_pts);
            sta_assert_error(next>=0);

            index=sphere_NN[current*num_sphere_pts+std::rand()%next]-1;
        }
        //printf("%d %d\n",index,num_sphere_pts);

        sta_assert_error(index>=0);
        if ((index>=num_sphere_pts)||(index<0))
            printf("%d %d\n",index,num_sphere_pts);
        sta_assert_error(index<num_sphere_pts);
// 	  n[2]=sphere_pts[3*index+2];
// 	  n[1]=sphere_pts[3*index+1];
// 	  n[0]=sphere_pts[3*index];


        n[0]=sphere_pts[3*index+2];
        n[1]=sphere_pts[3*index+1];
        n[2]=sphere_pts[3*index];


    }

    sta_assert_error(index!=-1);
    return index;
}








#endif
