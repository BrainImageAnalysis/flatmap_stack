/*#############################################################################
 *
 *	Copyright 2011 by Henrik Skibbe and Marco Reisert
 * 
 *     
 *	This file is part of the STA-ImageAnalysisToolbox
 * 
 *	STA-ImageAnalysisToolbox is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 * 
 *	STA-ImageAnalysisToolbox is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 * 
 *	You should have received a copy of the GNU General Public License
 *	along with STA-ImageAnalysisToolbox. 
 *	If not, see <http://www.gnu.org/licenses/>.
 *
 *
*#############################################################################*/

#ifndef MHS_COMMON_ERROR_H
#define MHS_COMMON_ERROR_H
#include <ctime>
#include <list>
#include <sstream>
#include <string>
#include <iomanip>
#include <assert.h>
#include <list>
#include <string>
#include <cstddef>
#include <complex>
#include <cmath>
#include <sstream>
#include <cstddef>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <limits>
#include <sys/stat.h>


#ifdef __linux__
  #include <sys/time.h>
#endif


//mxAssert sometimes didn't work, so we use this workaround

// #ifdef _SUPPORT_MATLAB_ 
//   #define sta_Assert( S ) \
//   { \
//     if (!( S )) \
//     {\
//      std::stringstream s; \
//      s<<"Assertion \""<<#S<<"\" failed (line "<<__LINE__<<" in "<<__FILE__<<")";\
//      printf("%s\n",s.str().c_str()); \
//       mexErrMsgTxt(s.str().c_str()); \
//     }\
//   }
// #else
//   #define sta_Assert( S ) \
//   { \
//     if (!( S )) \
//     {\
//       assert(S); \
//     }\
//   }
// #endif

#ifdef _SUPPORT_MATLAB_ 
  #define sta_assert( S ) \
  { \
    if (!( S )) \
    {\
     std::stringstream s; \
     s<<"Assertion \""<<#S<<"\" failed (line "<<__LINE__<<" in "<<__FILE__<<")";\
     printf("%s\n",s.str().c_str()); \
      mexErrMsgTxt(s.str().c_str()); \
    }\
  }
#else
  #define sta_assert( S ) \
  { \
    if (!( S )) \
    {\
      assert(S); \
    }\
  }
#endif


#define sta_assert_error( S ) \
{ \
  if (!( S )) \
  {\
    std::stringstream s; \
    s<<"Assertion \""<<#S<<"\" failed (line "<<__LINE__<<" in "<<__FILE__<<")";\
    printf("%s\n",s.str().c_str()); \
    throw mhs::STAError(s.str());\
  }\
}

#define sta_assert_error_m( S , M ) \
{ \
  if (!( S )) \
  {\
    std::stringstream s; \
    s<<"Assertion \""<<#S<<"\" failed (line "<<__LINE__<<" in "<<__FILE__<<")"<<std::endl<<"Message: "<<M;\
    printf("%s\n",s.str().c_str()); \
    throw mhs::STAError(s.str());\
  }\
}

#define sta_assert_error_c( S , C ) \
{ \
  if (!( S )) \
  {\
    std::stringstream s; \
    s<<"Assertion \""<<#S<<"\" failed (line "<<__LINE__<<" in "<<__FILE__<<")";\
    C;\
    printf("%s\n",s.str().c_str()); \
    throw mhs::STAError(s.str());\
  }\
}



//#define _DEBUG_CHECKS_0_
/*
 _DEBUG_CHECKS_0_ // fast, minor checks
 _DEBUG_CHECKS_1_ // further checks
 */


#ifdef _DEBUG_CHECKS_1_ 
  #define _DEBUG_CHECKS_0_
#endif

#ifdef _DEBUG_CHECKS_0_
  #define sta_assert_debug0( S ) \
  { \
    if (!( S )) \
    {\
      std::stringstream s; \
      s<<"Assertion \""<<#S<<"\" failed (line "<<__LINE__<<" in "<<__FILE__<<")";\
      printf("%s\n",s.str().c_str()); \
      throw mhs::STAError(s.str());\
    }\
  }
#else
  #define sta_assert_debug0( S ) {}
#endif

namespace  mhs
{

  



bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}  
  
  
  /// function return value 
enum STA_RESULT {
    STA_RESULT_SUCCESS=0,
    STA_RESULT_FAILED=-1,
    STA_RESULT_POINTLIMIT=-2,      
};  

/// the STA error class
class STAError
{
public:
    STAError() {}

    STAError(const std::string& message, STA_RESULT error_code=STA_RESULT_FAILED)
            :_message(message),_error_code(error_code)
    {}
    
    STAError(STA_RESULT error_code)
            :_error_code(error_code)
    {_message="";}    


    template<class T>
    STAError & operator<<(const T &data)
    {
        std::ostringstream os;
        os << data;
        _message += os.str();
        return *this;
    }
    
    void setCode(STA_RESULT error_code)
    {
      _error_code=error_code;
    }
    
    STA_RESULT getCode() const
    {
      return _error_code;
    }    

    /// \returns  the error string
    std::string str() const
    {
        return _message;
    }


    /// \returns  the error c-string
    const char* what() const
    {
        return _message.c_str();
    }

private:
    std::string _message;
    STA_RESULT _error_code;
};


class CtimeStopper
{
private:
    double _startSecond;
    double _lastSecond;
    double _lastSecond_request;
    int steps;
    std::string startevent;
    std::list<std::string> event_names;
    std::list<double> event_times;
    std::list<double> event_all_times;
public:
    CtimeStopper(std::string event="")
    {
        startevent=event;
        steps=0;
#ifdef __linux__	
        struct timeval _Time;
        gettimeofday(&_Time,NULL);
        _startSecond=_Time.tv_sec + 0.000001*_Time.tv_usec;
#else
	_startSecond=0;
#endif
        _lastSecond=_startSecond;
	_lastSecond_request=_startSecond;
	if (startevent!="")
	{
	  addEvent(startevent);
	  //printEvents();
	}	
    }
//     ~CtimeStopper()
//     {
//        if (startevent!="")
//        {
// 	 addEvent("done");
// 	 printEvents();
//        }
//     }
    
    
    double get_total_runtime()
    {
#ifdef __linux__      
        struct timeval  currentTime;
        gettimeofday(&currentTime,NULL);
        double currentSeconds=currentTime.tv_sec + 0.000001*currentTime.tv_usec;
#else
	double currentSeconds=0;
#endif
	return currentSeconds-_startSecond;
    }
    
     double get_last_call_runtime()
    {
#ifdef __linux__      
        struct timeval  currentTime;
        gettimeofday(&currentTime,NULL);
        double currentSeconds=currentTime.tv_sec + 0.000001*currentTime.tv_usec;
#else
	double currentSeconds=0;
#endif
	double time_diff=currentSeconds-_lastSecond_request;
	_lastSecond_request=currentSeconds;
	return time_diff;
    }
    
    void addEvent(std::string event)
    {
#ifdef __linux__      
        struct timeval  currentTime;
        gettimeofday(&currentTime,NULL);
        double currentSeconds=currentTime.tv_sec + 0.000001*currentTime.tv_usec;
#else
	double currentSeconds=0;
#endif
	
	event_times.push_back(currentSeconds-_lastSecond);
	event_all_times.push_back(currentSeconds-_startSecond);
	event_names.push_back(event);	
        _lastSecond=currentSeconds;
        steps++;
    }    
    
    void printEvents(int precision=3)
    {
      printf("\n computation time (summary):\n");
      printf("    all         step\n");
      std::list<std::string>::iterator iter2=event_names.begin();
      std::list<double>::iterator iter3=event_all_times.begin();
      for (std::list<double>::iterator iter=event_times.begin();
	  iter!=event_times.end();iter++,iter2++,iter3++)
      {
	std::stringstream s;
	s<<std::fixed<<std::setprecision(precision)<<*iter3<<" sec. | "<<*iter<<" sec. |\t("<<*iter2<<")";
	printf("%s\n",s.str().c_str());
      }
      printf("\n");
    }
};


}



#endif

