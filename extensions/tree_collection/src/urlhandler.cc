#include <sstream>
#include <iostream>
#include "urlhandler.h"

using namespace std;

int32_t UrlHandler::processURL(const char* url) {
   request_.SetURL(url);
   request_.SetMethod("GET");
   request_.SetFollowRedirects(true);

   download_complete_ = false;
   had_error_ = false;

   pp::CompletionCallback cc = factory_.NewCallback(&UrlHandler::DidCompleteIO);
   int32_t rv = loader_.Open(request_, cc);
   if (rv != PP_OK_COMPLETIONPENDING)
      cc.Run(rv);

   return rv;
}

void UrlHandler::DidCompleteIO(int32_t result) {
   if (result > 0) {
      // buf_ now contains 'result' number of bytes from the URL.
      ProcessBytes(buf_, result);
      ReadMore();
   } else if (result == PP_OK && !did_open_) {
      // Headers are available, and we can start reading the body.
      did_open_ = true;
      ProcessResponseInfo(loader_.GetResponseInfo());
      ReadMore();
   } else if (result == PP_OK) {
      // Done reading
      std::cerr << "Done reading URL" << std::endl;
      if(name_ != "") {
         instance_->PostMessage(pp::Var(name_+": OK"));
      }
      download_complete_ = true;
   } else {
      // Error
      std::cerr << "Error reading URL: " << result << std::endl;
      if(name_ != "") {
         instance_->PostMessage(pp::Var(name_+": FAILED"));
      }
      had_error_ = true;
   }
}

void UrlHandler::ReadMore() {
   pp::CompletionCallback cc = factory_.NewCallback(&UrlHandler::DidCompleteIO);
   int32_t rv = loader_.ReadResponseBody(buf_, sizeof(buf_), cc);
   if (rv != PP_OK_COMPLETIONPENDING)
      cc.Run(rv);
}

void UrlHandler::ProcessResponseInfo(const pp::URLResponseInfo& response_info) {
   std::cerr << response_info.GetStatusCode() << std::endl;
   std::cerr << response_info.GetStatusLine().AsString() << std::endl;
}

void UrlHandler::ProcessBytes(const char* bytes, int32_t length) {
   ss_.write(bytes,length);
}

std::string UrlHandler::getData() {
   return ss_.str();
}

