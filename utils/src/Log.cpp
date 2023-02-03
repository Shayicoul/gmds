/*----------------------------------------------------------------------------*/
/** \file    Log.cpp
 *  \author  F. LEDOUX
 *  \date    09/17/2009
 */
/*----------------------------------------------------------------------------*/
#include "../inc/gmds/utils/Log.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
LogLevel Log::reporting_level_ = LOG_INFO;
/*----------------------------------------------------------------------------*/
Log* Log::mng_=0;
/*----------------------------------------------------------------------------*/
Log::Log(){}
/*----------------------------------------------------------------------------*/
Log::~Log(){}
/*----------------------------------------------------------------------------*/
LogLevel& Log::reportingLevel(){return reporting_level_;}
/*----------------------------------------------------------------------------*/
void Log::addStream(LogStream& AStream) {out_.push_back(AStream);}
/*----------------------------------------------------------------------------*/
Log& Log::mng(){
  if(mng_)
    return *mng_;

  mng_ = new Log();
  return *mng_;
}
/*----------------------------------------------------------------------------*/
