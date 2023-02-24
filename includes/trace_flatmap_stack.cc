#include "mex.h"

#include <unistd.h>
#include <complex>
#include <map>
#include "mhs_error.h"
#include "mhs_vector.h"
#include <algorithm> 


#define _SUPPORT_MATLAB_
#include "sta_mex_helpfunc.h"
// #define _SUPPORT_MATLAB_ 

#include <complex>
#include <vector>
#include <string>
#include <ctime>
#include <list>
#include <sstream>
#include <string>
#include <limits>
#include <omp.h>

// #define SUB2IND(X, Y, Z, shape)  (((Z)*(shape[1])+(Y)) *(shape[2])+(X)) 
#define SUB2IND(X, Y, Z, shape)  (((Z)*(shape[1])+(Y)) *(shape[0])+(X)) 
#include "mhs_vector.h"


template<typename T>
void hsv2rgb ( T h, T s ,T v,T &r, T & g ,T & b )
{
    T      hh, p, q, t, ff;
    int        i;


    if ( s <= 0.0 ) {  
        r = v;
        g = v;
        b = v;
    }
    hh = h;
    if ( hh >= 360.0 ) {
        hh = 0.0;
    }
    hh /= 60.0;
    i = ( int ) hh;
    ff = hh - i;
    p = v * ( 1.0 - s );
    q = v * ( 1.0 - ( s * ff ) );
    t = v * ( 1.0 - ( s * ( 1.0 - ff ) ) );

    switch ( i ) {
    case 0:
        r = v;
        g = t;
        b = p;
        break;
    case 1:
        r = q;
        g = v;
        b = p;
        break;
    case 2:
        r = p;
        g = v;
        b = t;
        break;

    case 3:
        r = p;
        g = q;
        b = v;
        break;
    case 4:
        r = t;
        g = p;
        b = v;
        break;
//     case 5:
    default:
        r = v;
        g = p;
        b = q;
        break;
    }
}



bool ind2sub(std::size_t * shape, std::size_t indx, std::size_t &X, std::size_t &Y,std::size_t &Z)
{
    if (indx>=shape[0]*shape[2]*shape[1])
        return false;
    
Z = indx / (shape[0]*shape[1]);
     indx %=  shape[0]*shape[1];
//     
     Y = indx / (shape[0]);
     X = indx % shape[0];
//       Z = indx / (shape[2]*shape[1]);
//      indx %=  shape[2]*shape[1];
// //     
//      Y = indx / (shape[1]);
//      X = indx % shape[1];
    
    
     return true;   
}

//mex trace_flatmap_stack.cc  -largeArrayDims -lgomp CXXFLAGS=" -O3   -Wfatal-errors  -std=c++11 -fopenmp-simd -fopenmp  -fPIC -march=native"



template <typename T>
void _mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  try
        {
  
              
    const mxArray * params=prhs[nrhs-1];

    mxArray * ofield = NULL;
    mxArray * mask = NULL;
    mxArray * seeds = NULL;    
    mxArray * dist = NULL; 

    int mode = 3;
    std::size_t max_search = 200;
    
    std::size_t rad = 1;
    T density = 3;
    T stepwidth = 1;



    if ( mxIsCell ( params  )  )
    {
        if ( mhs::mex_hasParam ( params,"max_search" ) !=-1 )
        {
            max_search=mhs::mex_getParam<std::size_t > ( params,"max_search",1 ) [0];
        }
        
         if ( mhs::mex_hasParam ( params,"stepwidth" ) !=-1 )
        {
            stepwidth=mhs::mex_getParam<T> ( params,"stepwidth",1 ) [0];
        }
        
                 if ( mhs::mex_hasParam ( params,"density" ) !=-1 )
        {
            density=mhs::mex_getParam<T> ( params,"density",1 ) [0];
        }
        
        if ( mhs::mex_hasParam ( params,"rad" ) !=-1 )
        {
            rad=mhs::mex_getParam<std::size_t > ( params,"rad",1 ) [0];
        }
        
        if ( mhs::mex_hasParam ( params,"mode" ) !=-1 )
        {
            mode=mhs::mex_getParam<int> ( params,"mode",1 ) [0];
        }
        
        if ( mhs::mex_hasParam ( params,"ofield" ) !=-1 )
        {
            ofield=mhs::mex_getParamPtr( params,"ofield");
        }
        
        if ( mhs::mex_hasParam ( params,"mask" ) !=-1 )
        {
            mask=mhs::mex_getParamPtr( params,"mask");
        }

        if ( mhs::mex_hasParam ( params,"dist" ) !=-1 )
        {
            dist=mhs::mex_getParamPtr( params,"dist");
        }
        
         if ( mhs::mex_hasParam ( params,"seeds" ) !=-1 )
        {
            seeds=mhs::mex_getParamPtr( params,"seeds");
        }
    }

    sta_assert_error(ofield!=NULL);
    sta_assert_error(mask!=NULL);
    sta_assert_error(seeds!=NULL);
    sta_assert_error(dist!=NULL);

    mhs::dataArray<T>  ofield_;
    ofield_=mhs::dataArray<T> (ofield);	 
    sta_assert_error(ofield_.dim.size()==4);
    printf("%d %d %d %d\n",ofield_.dim[0],ofield_.dim[1],ofield_.dim[2],ofield_.dim[3]);
    sta_assert_error(ofield_.dim[3]==3);        
    sta_assert_error(ofield_.dim[0]>1);
    T * ofield_p = ofield_.data;
    
    
    mhs::dataArray<T>  mask_;
    mask_ = mhs::dataArray<T> (mask);	 
    sta_assert_error(mask_.dim.size()==3);    
    sta_assert_error(mask_.dim[0]>1);
    T * mask_p = mask_.data;    

    mhs::dataArray<T>  dist_;
    dist_ = mhs::dataArray<T> (dist);	 
    sta_assert_error(dist_.dim.size()==3);    
    sta_assert_error(dist_.dim[0]>1);
    T * dist_p = dist_.data;       
    
    mhs::dataArray<T>  seeds_;
    seeds_ = mhs::dataArray<T> (seeds);	 
    sta_assert_error(seeds_.dim.size()==2);    
    sta_assert_error(seeds_.dim[0]>1);
    T * seeds_p = seeds_.data;  

    std::size_t shape[3];
    shape[0] = mask_.dim[2];
    shape[1] = mask_.dim[1];
    shape[2] = mask_.dim[0];    
    printf("shape ; %d %d %d\n",shape[0],shape[1],shape[2]);    
    //   mhs::dataArray<T>  gene_p_;
//      T * gene_p = NULL;
  
    
    
    printf("%d %d %d %d\n",ofield_.dim[0],ofield_.dim[1],ofield_.dim[2],ofield_.dim[3]);
    printf("%d %d %d\n",mask_.dim[0],mask_.dim[1],mask_.dim[2]);
    printf("%d %d\n",seeds_.dim[0],seeds_.dim[1]);
    printf("%d %d %d\n",dist_.dim[0],dist_.dim[1],dist_.dim[2]);

    T *result = NULL; 
    /*
    mwSize ndims[3];
    ndims[2]=mask_.dim[0];
    ndims[1]=mask_.dim[1];
    ndims[0]=mask_.dim[2];
    plhs[0] = mxCreateNumericArray ( 3,ndims,mhs::mex_getClassId<T>(),mxREAL );
    */
    mwSize ndims[4];
    ndims[3]=mask_.dim[0];
    ndims[2]=mask_.dim[1];
    ndims[1]=mask_.dim[2];
    ndims[0]=3;
    plhs[0] = mxCreateNumericArray ( 4,ndims,mhs::mex_getClassId<T>(),mxREAL );
    result = ( T * ) mxGetData ( plhs[0]);;

    std::size_t nvox = shape[0]*shape[1]*shape[2];
    std::size_t seed_stride = seeds_.dim[1];

    plhs[1] = mxCreateCellMatrix(seed_stride, 1);


    
    
    for (std::size_t t =0; t<seeds_.dim[1];t++)
    {
        //bool drawit = myrand(100.0)>90;
        bool drawit = myrand(100.0)>(100-density);
        
        //if (myrand(100.0)>3)
        //    continue;
        T r;
        T g;
        T b;
        //hsv2rgb ( T h, T s ,T v,T &r, T & g ,T & b )
        hsv2rgb(myrand(360.0),1.0f,1.0f,r,g,b);
        
        std::list<Vector<float,4> > trackt;
        //int signi = 1;
        for (int signi=0;signi<2;signi++)
        {
            Vector<T,3> seed(seeds_p[t],seeds_p[t+seed_stride],seeds_p[t+2*seed_stride]);

            float sign = signi>0 ? -1 : 1;
            //seed.print();
            int count = max_search;
            bool ok = true;
            bool init = false;
            while (count>0 && ok)
            {
                count--;
                int x = std::round(seed[0]);
                int y = std::round(seed[1]);
                int z = std::round(seed[2]);

                std::size_t pIND = SUB2IND(x,y,z,shape);
                                        
                if (pIND>0 && pIND<nvox)
                {
                    if (mask_p[pIND]>0.5)
                    {
                       // Vector<float,3> pt;
                       // pt[0] = seed[0];
                       // pt[1] = seed[1];
                        //pt[2] = seed[2];
                        
                        Vector<float,4> step_;
                        step_[0] = seed[0];
                        step_[1] = seed[1];
                        step_[2] = seed[2];
                        step_[3] = dist_p[pIND];
                        if (signi==0)
                        {
                            trackt.push_back(step_);
                        }
                        else
                            if (init)
                            {
                                trackt.push_front(step_);
                            }
                        Vector<T,3> vec(ofield_p[3*pIND],ofield_p[3*pIND+1],ofield_p[3*pIND+2]);
                        //printf("%f %f %f\n",vec[0],vec[1],vec[2]);
                        vec.normalize();
                        //printf("%f %f %f\n",vec[0],vec[1],vec[2]);
                        seed += vec*(sign*0.5);

                        if (drawit)
                        {
                            result[3*pIND] = r;
                            result[3*pIND+1] = g;
                            result[3*pIND+2] = b;
                        }
                        //result[pIND*2] = 0.5;
                        //result[pIND*3] = 0.1;
                    }else
                    {
                        ok = false;
                    }
                }else
                {
                    printf("oo image\n");
                };
                init = true;
            }
            
        }
        //printf("yea %d\n",seed_stride);
        //break;
       
         mwSize ndims[2];
         ndims[0]=trackt.size();
         ndims[1]=4;
         mxArray * trc = mxCreateNumericArray ( 2,ndims,mhs::mex_getClassId<T>(),mxREAL );
         T * tr = ( T * ) mxGetData ( trc);
         mxSetCell(plhs[1],t,trc);
         std::size_t c = 0;
         int t_size = trackt.size();
         for (std::list<Vector<float,4> >::iterator iter = trackt.begin();iter!=trackt.end();iter++)
         {
                Vector<float,4> & p = * iter;
                tr[c] = p[0];
                tr[c+t_size] = p[1];
                tr[c+2*t_size] = p[2];
                tr[c+3*t_size] = p[3];
                c++;
         }
         
    }

   

    /*
    for (int i=0;i<3;i++)
    {
        sta_assert_error(cmask_s_.dim[i]==bmask_o_.dim[i]);
        sta_assert_error(bmask_n_.dim[i]==bmask_o_.dim[i]);
    }
    */

    
        
     /*
     while (count>0 && ok<2)
        {
            count--;
            int x = std::round(seed[0]);
            int y = std::round(seed[1]);
            int z = std::round(seed[2]);

            std::size_t pIND = SUB2IND(x,y,z,shape);
                                    
            if (pIND>0 && pIND<nvox)
            {

                Vector<T,3> vec(ofield_p[3*pIND],ofield_p[3*pIND+1],ofield_p[3*pIND+2]);
                vec.normalize();
                seed += vec*0.5;

                if (mask_p[pIND]>0.5)
                {
                    if (ok==0)
                    {
                        ok = 1;
                    }

                    result[3*pIND] = r;
                    result[3*pIND+1] = g;
                    result[3*pIND+2] = b;
                    //result[pIND*2] = 0.5;
                    //result[pIND*3] = 0.1;
                }else
                if (ok>0)
                {
                    ok++;
                }
            }else
            {
                printf("oo image\n");
            };
        }*/
        
        



    }
    catch ( mhs::STAError & error )
    {
        printf ( "error cleaning up!!\n" );
        mexErrMsgTxt ( error.what() );
    }
    catch ( ... )
    {
        printf ( "error cleaning up!!\n" );
        mexErrMsgTxt ( "exeption ?" );
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{    
 

    if (nlhs<2)
        mexErrMsgTxt ( "error: plhs<2\n" );

    if ( nrhs<1 )
        mexErrMsgTxt ( "error: nrhs<1\n" );
     _mexFunction<float> ( nlhs, plhs,  nrhs, prhs );
   
}

