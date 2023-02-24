#ifndef STA_MEX_HELPFUNC
#define STA_MEX_HELPFUNC

#include "mex.h"
#include "mhs_error.h"
#include <complex>
#include <vector>
#include <algorithm> 
#include <string>
#include "omp.h"
#include <stdint.h> 



#ifndef printf
  #define printf mexPrintf
#endif

  #include <stdio.h>  /* defines FILENAME_MAX */
  
  
// #ifdef _WIN32
//   #include <stdint.h>
//   typedef uint8_t u_int8_t;
//   typedef uint16_t u_int16_t;
//   typedef uint32_t u_int32_t;
//   typedef uint64_t u_int64_t;
// #else
//   #ifdef __APPLE__
// 
//   #endif
// #endif
  
// #ifdef WINDOWS
//     #include <direct.h>
//     #define GetCurrentDir _getcwd
// #else
//     #include <unistd.h>
//     #define GetCurrentDir getcwd
//  #endif


// #ifdef __cplusplus 
//     extern "C" bool utIsInterruptPending();
// #else
//     extern bool utIsInterruptPending();
// #endif

//  #ifndef _mxGetPropertyPtr_
//   #define _mxGetPropertyPtr_
//     /*
//      * mxGetPropertyPtr by James Tursa
//      * see 
//      * mxGetPropertyPtr/mxGetPropertyPtr_20110307.pdf
//      * for further details
//     */
//     #include "./mxGetPropertyPtr/mxGetPropertyPtr.c"
//   #endif 




namespace  mhs
{

// /*
//  char cCurrentPath[FILENAME_MAX];
// 
//  if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
//   {
//   return errno;
//   }
// 
// cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */
// 
// printf ("The current working directory is %s", cCurrentPath);*/
//   

template <typename T>
std::string mex_getClassString()
{
    return  "unknown";
}

template<>
std::string mex_getClassString<float>()
{
    return  "float";
}

template<>
std::string mex_getClassString<double>()
{
    return  "double";
}

template<>
std::string mex_getClassString<bool>()
{
    return  "bool";
}
  
template<>
std::string mex_getClassString<uint8_t>()
{
    return  "uint8_t";
}

template<>
std::string mex_getClassString<uint32_t>()
{
    return  "uint32_t";
}
   
template<>
std::string mex_getClassString<uint64_t>()
{
    return  "uint64_t";
}
   
  
template <typename T>
mxClassID mex_getClassId()
{
    return  mxUNKNOWN_CLASS;
}

template<>
mxClassID mex_getClassId<bool>()
{
    return mxLOGICAL_CLASS;
}


template<>
mxClassID mex_getClassId<float>()
{
    return mxSINGLE_CLASS;
}

template<>
mxClassID mex_getClassId<double>()
{
    return  mxDOUBLE_CLASS;
}


template<>
mxClassID mex_getClassId<uint8_t>()
{
    return   mxUINT8_CLASS;
}


template<>
mxClassID mex_getClassId<uint64_t>()
{
    return   mxUINT64_CLASS;
}

template<>
mxClassID mex_getClassId<uint32_t>()
{
    return   mxUINT32_CLASS;
}  
 
template <typename T>
bool mex_israwarray(const mxArray * array)
{
  return false;
}

template <>
bool mex_israwarray<float>(const mxArray * array)
{
  return mxIsClass(array,"single");
}

template <>
bool mex_israwarray<double>(const mxArray * array)
{
  return mxIsClass(array,"double");
}

template <>
bool mex_israwarray<uint8_t>(const mxArray * array)
{
  return mxIsClass(array,"uint8");
}


template <typename T>
mxArray * mex_create_value(T value)
{
      mwSize ndim=1;
      mxArray *handle = mxCreateNumericArray(1,&ndim,mex_getClassId<T>(),mxREAL);
     T * handle_p=(T  *)mxGetPr(handle);
     *handle_p=value;
     return handle;
}


 
template<typename T>
T * mex_getPtr(const mxArray *mx_array,
	       std::vector<std::size_t> & dim)
{   
    //const mxArray *mx_array= mxGetField(feature_struct,0,(char *)(fieldname.c_str()));
    sta_assert_error(mx_array!=NULL);
    if (!mex_israwarray<T>(mx_array))
    {
      mhs::STAError error;
      error << "expecting data type " <<mex_getClassString<T>()<<", but it is "<<mxGetClassName( ( mx_array ))<<"\n";
      throw error;
    }
    
    
    const std::size_t  mx_ndims = mxGetNumberOfDimensions ( mx_array );
    const mwSize *mx_dims = mxGetDimensions ( mx_array );
    dim.resize(mx_ndims);
    for (std::size_t i=0;i<mx_ndims;i++)
      dim[i]=mx_dims[i];
    if (mxGetClassID ( mx_array ) !=mex_getClassId<T>())
    {
      mhs::STAError error;
      error << "data type" <<mex_getClassId<T>()<< " expected, but it is "<<mxGetClassName( ( mx_array ))<<"\n";
      //error << "wrong data type\n";
      throw error;
    }
    
    return  ( T * ) ( mxGetData ( mx_array ) );
}

template<typename T>
T * mex_getPtr(const mxArray *mx_array)
{   
   try
   {
      std::vector<std::size_t>  dim;
      return  mex_getPtr<T>(mx_array,dim);
    }catch (mhs::STAError error)
    {
	throw error;
    }
}


template<typename T>
T * mex_getFieldPtr(const mxArray *mx_array,
	       std::string fieldname,
	       std::vector<std::size_t> & dim,int indx=0)
{   
    sta_assert_error(mxIsStruct(mx_array));
    const mxArray *mx_array2= mxGetField(mx_array,indx,(char *)(fieldname.c_str()));
    T * ptr=NULL;
    try
    {
      ptr=mex_getPtr<T>(mx_array2,dim);
    }catch (mhs::STAError error)
    {
	error << " field : "<<fieldname<<"\n";
	throw error;
    }
     
    return ptr;
}

template<typename T>
T * mex_getFieldPtr(const mxArray *mx_array,
	       std::string fieldname)
{   
    std::vector<std::size_t>  dim;
    T * ptr=NULL;
    try
    {
      ptr=mex_getFieldPtr<T>(mx_array,fieldname,dim);
    }catch (mhs::STAError error)
    {
	throw error;
    }
    return ptr;
}

template<typename T>
T * mex_getFieldPtrCreate(mxArray *mx_array,mwSize dim,std::string fieldname,mwSize indx=0)
{
  sta_assert_error(mxIsStruct(mx_array));
  sta_assert_error(indx<mxGetNumberOfElements(mx_array));
  T * ptr=NULL;
  
  bool field_exists=mxGetFieldNumber(mx_array,(char *)(fieldname.c_str()))>-1;
  
  mxArray *struct_array=NULL;
   try
    {
      if (!field_exists)
      {
// 	printf("E0\n");
	sta_assert_error(mxAddField(mx_array, (char *)(fieldname.c_str()))>-1);
      }
      sta_assert_error(mxGetFieldNumber(mx_array,(char *)(fieldname.c_str()))>-1);
//       printf("E1\n");
      struct_array=mxGetField(mx_array,indx,(char *)(fieldname.c_str()));
      if (struct_array==NULL)
      {
// 	printf("E2\n");
	//struct_array= mxCreateNumericArray ( 1,&dim,mhs::mex_getClassId<T>(),mxREAL);
 	
        mxArray * pdata = mxCreateNumericArray ( 1,&dim,mhs::mex_getClassId<T>(),mxREAL);
// 	printf("E3\n");
 	mxSetField(mx_array,indx,(char *)(fieldname.c_str()),pdata);  
	struct_array=mxGetField(mx_array,indx,(char *)(fieldname.c_str()));
      }
//       printf("E4\n");
      sta_assert_error(struct_array!=NULL);
      sta_assert_error(indx<mxGetNumberOfElements(mx_array));
      sta_assert_error(dim==mxGetNumberOfElements(struct_array));
    
//       printf("E5\n");
      ptr=mex_getPtr<T>(struct_array);
      
  }catch (mhs::STAError error)
  {
    throw error; 
  }
  printf("E6\n");
    return ptr;
}


void copyField(const mxArray * from,mxArray * to,std::string fieldname,int indxFrom=0, int indxTo=0 )
{
	  sta_assert_error(mxIsStruct(from));
	  sta_assert_error(mxIsStruct(to));
	  sta_assert_error(mxGetField(from,indxFrom,(char *)(fieldname.c_str()))!=NULL);
	  sta_assert_error(mxGetField(to,indxTo,(char *)(fieldname.c_str()))==NULL);
	  mxAddField(to,(char *)(fieldname.c_str()));
	  mxSetField(to, indxTo,(char *)(fieldname.c_str()), mxDuplicateArray(mxGetField(from,indxFrom,(char *)(fieldname.c_str()))));
}

void copyToField(const mxArray * from,mxArray * to,std::string fieldname, int indxTo=0 )
{
	  sta_assert_error(mxIsStruct(to));
	  sta_assert_error(mxGetField(to,indxTo,(char *)(fieldname.c_str()))==NULL);
	  mxAddField(to,(char *)(fieldname.c_str()));
	  mxSetField(to, indxTo,(char *)(fieldname.c_str()), mxDuplicateArray(from));
}

mxArray * createField(mxArray * mx_array,std::string fieldname,int indx=0)
{
    sta_assert_error(mxAddField(mx_array, (char *)(fieldname.c_str()))>-1);
    sta_assert_error(mxGetFieldNumber(mx_array,(char *)(fieldname.c_str()))>-1);
    return mxGetField(mx_array,indx,(char *)(fieldname.c_str()));
}
      
      

    


template<typename T>
class dataArray
{
public:
  T * data;
  std::vector<std::size_t>  dim;
  dataArray(const mxArray *mx_array,std::string field="",int indx=0)
  {
    try {
      if (mxIsStruct(mx_array))
      {
	data=mex_getFieldPtr<T>(mx_array,field,dim,indx);
      }else
      {
	data=mex_getPtr<T>(mx_array,dim);
      }
      std::reverse(dim.begin(),dim.end());
    }catch (mhs::STAError error)
    {
	throw error;
    }
  }
  
  dataArray() : data(NULL){};
  
  
  std::size_t get_num_elements()
  {
    std::size_t  elements=1;
    for (std::size_t i=0;i<dim.size();i++)
    {
      elements*=dim[i];
    }
    if (dim.size()==0)
      return 0;
    return elements;
    
  }
  
  void print()
  {
    printf("[");
    for (int i=0;i<dim.size();i++)
    printf(" %d ",dim[i]);  
    printf("]\n");
  }
};  
  

void mex_initMatlabFFTW()
{
    mexEvalString("ifft(fft(single([1,2,3])));ifft(fft(double([1,2,3])));");
}

void mex_dumpStringNOW()
{
#ifndef _VALGIND_CONSOLE_  
    int thread_id=omp_get_thread_num();
    if (thread_id==0)
    {
    mexEvalString("drawnow;");
    }
#endif
}




template <typename T>
std::complex<T> mex_getComplexValue(const mxArray *  ptr)
{
    if (mxGetNumberOfDimensions(ptr)==2)
    {
        if (mxGetDimensions(ptr)[0]*mxGetDimensions(ptr)[1]==2)
            return std::complex<T>((*(T*)mxGetData(ptr)),(*((T*)mxGetData(ptr)+1)));
        if (mxGetDimensions(ptr)[0]*mxGetDimensions(ptr)[1]==1)
            return (*(T*)mxGetData(ptr));
    }

    mexErrMsgTxt("cannot parse the alpha parameter?!\n");
    return T(1.0);
}


std::string mex_mex2string(const mxArray * s)
{
     char *  S=mxArrayToString(s);
     if (S!=NULL)
         return(std::string(S));
     else mexErrMsgTxt("error parsing string\n");
     return("");


//     char buffer[255];
//     if (mxGetString(s, buffer, 255)==0)
//         return(std::string(buffer));
//     else mexErrMsgTxt("error parsing string\n");
//     return("");
}


mxArray * mex_string2mex(std::string s)
{
  return mxCreateString(s.c_str());
}

template<typename T>
std::string mex_value2string(T value)
{
    std::stringstream s;
    s<<value;
    return s.str();
}



template <typename T>
std::string mex_getPrecisionFlag()
{
    return  "unknown";
}

template<>
std::string mex_getPrecisionFlag<float>()
{
    return "single";
}

template<>
std::string mex_getPrecisionFlag<double>()
{
    return  "double";
}



template<typename T>
bool mex_isStaFieldStruct(const mxArray * s)
{
    if (!mxIsStruct(s))
        return false;
    if (mxGetField(s,0,"storage")==NULL)
        return false;
    if (mxGetField(s,0,"L")==NULL)
        return false;
    if (mxGetField(s,0,"type")==NULL)
        return false;
    if (mxGetField(s,0,"data")==NULL)
        return false;
    if (mxGetClassID(mxGetField(s,0,(char*)"data"))!=mex_getClassId<T>())
        return false;
    return true;
}



template<typename T>
std::vector<T> mex_getParam_helper(std::string params)
{
    for (std::size_t a=0;a<params.length();a++)
        if (params[a]==',')
            params[a]=' ';
    std::vector<T> param_v;
    std::stringstream param_s;
    param_s<<params;
    int eexit=0;
    while (!param_s.eof() && param_s.good())
    {
        T tmp;
        param_s>>tmp;
        param_v.push_back(tmp);
        eexit++;
        if (eexit>42)
            mexErrMsgTxt("ahhh, something went completely wrong while parsing parameters! ");
    }
    return param_v;
}


int mex_hasParam(const mxArray * plist,std::string s,bool required=false)
{
printf("mex_hasParam ???\n");    
    if (!mxIsCell(plist))
        mexErrMsgTxt("error hasParam: parameterlist is not a cell array\n");
printf("mex_hasParam ???\n");      
    if (mxGetNumberOfElements(plist)%2==1)
        mexErrMsgTxt("error hasParam:  length of parameterlist is not even");
printf("mex_hasParam 0\n");
    for (std::size_t t=0;t<mxGetNumberOfElements(plist);t+=2)
    {
        if (mxGetClassID(mxGetCell(plist,t))!=mxCHAR_CLASS)
        {
            std::stringstream s;
            s<<"error: parameter "<<t<<"is not a string"<<"\n";
            mexErrMsgTxt(s.str().c_str());
        }

        std::string param(mex_mex2string(mxGetCell(plist,t)));

	//printf("%s %s %d\n",s.c_str(),param.c_str(),(param==s));
	
        if (param==s) return (int)t;
    }
printf("mex_hasParam 1\n");
    if (required)
    {
        std::string error="error: missing parameter "+s+"\n";
        mexErrMsgTxt(error.c_str());
    }
printf("mex_hasParam 2\n");   
    return -1;
}



mxArray * mex_getParamPtr(const mxArray * plist,std::string s)
{
    int pos=mex_hasParam(plist,s);
    if (pos==-1)
        mexErrMsgTxt("error getParamPtr: couldn't find parameter\n");
    return(mxGetCell(plist,pos+1));
}



std::string  mex_getParamStr(const mxArray * plist,std::string s)
{
    if (mex_hasParam(plist,s)==-1)
        mexErrMsgTxt("error getParam: couldn't find parameter\n");

    mxArray * params=mex_getParamPtr(plist,s);



    if (mxGetClassID(params)!=mxCHAR_CLASS)
        mexErrMsgTxt("string parameter expected!\n");

    return(mex_mex2string(params));
}



template<typename T>
std::vector<std::complex<T> > mex_getParamC(const mxArray * plist,std::string s,int expected=-1)
{
    if (mex_hasParam(plist,s)==-1)
        mexErrMsgTxt("error getParam: couldn't find parameter\n");

    mxArray * params=mex_getParamPtr(plist,s);

    int numParams=mxGetNumberOfElements(params);
    if ((numParams!=expected)&&(expected!=-1))
        mexErrMsgTxt("number of parameters differs from number of expected parameters!\n");

    if (mxGetClassID(params)==mxCHAR_CLASS)
        mexErrMsgTxt("string parameter cannot be parsed with this function!\n");

    std::vector<std::complex<T> > result(numParams);


    for (int t=0;t<numParams;t++)
    {
        switch (mxGetClassID(params))
        {

        case mxSINGLE_CLASS:
        {
            float * data_r=(float *)mxGetPr(params);
            float * data_i=(float *)mxGetPi(params);
            if (data_i!=NULL)
                result[t]=std::complex<T>((T)data_r[t],(T)data_i[t]);
            else
                result[t]=(T)data_r[t];
        }
        break;

        case mxDOUBLE_CLASS:
        {
            double * data_r=(double *)mxGetPr(params);
            double * data_i=(double *)mxGetPi(params);
            if (data_i!=NULL)
                result[t]=std::complex<T>((T)data_r[t],(T)data_i[t]);
            else
                result[t]=(T)data_r[t];
        }
        break;
        default:
            mexErrMsgTxt("unsupported data type!\n");
        }
    }
    return result;
}

template<typename T>
std::vector<T> mex_getParam(const mxArray * plist,std::string s,int expected=-1)
{
    if (mex_hasParam(plist,s)==-1)
        mexErrMsgTxt("error getParam: couldn't find parameter\n");

    mxArray * params=mex_getParamPtr(plist,s);

    int numParams=mxGetNumberOfElements(params);
    if ((numParams!=expected)&&(expected!=-1))
        mexErrMsgTxt("number of parameters differs from number of expected parameters!\n");

    if (mxGetClassID(params)==mxCHAR_CLASS)
        mexErrMsgTxt("string parameter cannot be parsed with this function!\n");

    std::vector<T > result(numParams);



    for (int t=0;t<numParams;t++)
    {
        switch (mxGetClassID(params))
        {

        case mxSINGLE_CLASS:
        {
            float * data_r=(float *)mxGetPr(params);
            result[t]=(T)data_r[t];
        } break;

        case mxDOUBLE_CLASS:
        {
            double * data_r=(double *)mxGetPr(params);
            result[t]=(T)data_r[t];
        } break;

        case mxLOGICAL_CLASS:
        {
            bool * data_r=(bool *)mxGetPr(params);
            result[t]=(T)data_r[t];
        } break;
	 case  mxUINT64_CLASS:
        {
            uint64_t * data_r=(uint64_t *)mxGetPr(params);
            result[t]=(T)data_r[t];
        } break;
        default:
            mexErrMsgTxt("unsupported data type!\n");
        }
    }
    return result;
}



template<typename T>
T mex_get1Param(const mxArray * param)
{
    int numParams=mxGetNumberOfElements(param);
    if (numParams!=1)
        mexErrMsgTxt("number of parameters differs from 1 (number of expected parameters)!\n");
    T result=-1;
    switch (mxGetClassID(param))
    {

    case mxSINGLE_CLASS:
    {
        float * data_r=(float *)mxGetPr(param);
        result=(T)data_r[0];
    } break;

    case mxDOUBLE_CLASS:
    {
        double * data_r=(double *)mxGetPr(param);
        result=(T)data_r[0];
    } break;

    case mxLOGICAL_CLASS:
    {
        bool * data_r=(bool *)mxGetPr(param);
        result=(T)data_r[0];
    } break;

    default:
        mexErrMsgTxt("unsupported data type!\n");
    }
    return result;
}

template<typename T>
std::complex<T> mex_get1ParamC(const mxArray * param)
{
    int numParams=mxGetNumberOfElements(param);
    if (numParams!=1)
        mexErrMsgTxt("number of parameters differs from 1 (number of expected parameters)!\n");
    std::complex<T>  result=-1;
    switch (mxGetClassID(param))
    {

    case mxSINGLE_CLASS:
    {
        float * data_r=(float *)mxGetPr(param);
        float * data_i=(float *)mxGetPi(param);
        if (data_i!=NULL)
            result=std::complex<T>((T)data_r[0],(T)data_i[0]);
        else
            result=(T)data_r[0];
    }
    break;

    case mxDOUBLE_CLASS:
    {
        double * data_r=(double *)mxGetPr(param);
        double * data_i=(double *)mxGetPi(param);
        if (data_i!=NULL)
            result=std::complex<T>((T)data_r[0],(T)data_i[0]);
        else
            result=(T)data_r[0];
    }
    break;

    default:
        mexErrMsgTxt("unsupported data type!\n");
    }
    return result;
}
}

#endif





