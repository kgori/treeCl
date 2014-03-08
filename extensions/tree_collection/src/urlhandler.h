#include <string>
#include <sstream>

class UrlHandler {
   public:
      UrlHandler(pp::Instance* instance, const std::string name = "")
         : factory_(this),
         loader_(instance),
         request_(instance),
         did_open_(false),
         download_complete_(false),
         had_error_(false),
         name_(name),
         instance_(instance)
   {
   }
      int32_t processURL(const char* url);
      std::string getData();
      bool isFinished() { return download_complete_ || had_error_; }
      bool isComplete() { return download_complete_; }
      bool hadError() { return had_error_; }
   private:
      void DidCompleteIO(int32_t result);
      void ReadMore();
      void ProcessResponseInfo(const pp::URLResponseInfo& response_info);
      void ProcessBytes(const char* bytes, int32_t length);
      pp::CompletionCallbackFactory<UrlHandler> factory_;
      pp::URLLoader loader_;
      pp::URLRequestInfo request_;
      char buf_[4096];
      bool did_open_;
      bool download_complete_;
      bool had_error_;
      std::string name_;
      std::stringstream ss_;
      pp::Instance *instance_;
};

