/// @file darwintree.cc
/// This example demonstrates loading, running and scripting a very simple NaCl
/// module.  To load the NaCl module, the browser first looks for the
/// CreateModule() factory method (at the end of this file).  It calls
/// CreateModule() once to load the module code from your .nexe.  After the
/// .nexe code is loaded, CreateModule() is not called again.
///
/// Once the .nexe code is loaded, the browser than calls the CreateInstance()
/// method on the object returned by CreateModule().  It calls CreateInstance()
/// each time it encounters an <embed> tag that references your NaCl module.
///
/// The browser can talk to your NaCl module via the postMessage() Javascript
/// function.  When you call postMessage() on your NaCl module from the browser,
/// this becomes a call to the HandleMessage() method of your pp::Instance
/// subclass.  You can send messages back to the browser by calling the
/// PostMessage() method on your pp::Instance.  Note that these two methods
/// (postMessage() in Javascript and PostMessage() in C++) are asynchronous.
/// This means they return immediately - there is no waiting for the message
/// to be handled.  This has implications in your program design, particularly
/// when mutating property values that are exposed to both the browser and the
/// NaCl module.

#include <cstdio>
#include <string>

#include <cstdio>
#include <string>
#include <exception>
#include <iostream>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include <sstream>

#include "ProblemParser.h"
#include "MinSqTree.h"
#include "urlhandler.h"

namespace {
   // A method consists of a const char* for the method ID and the method's
   // declaration and implementation.
   const char * const kLoadMatricesMethodId = "loadMatrices";
   const char * const kLoadLabelsMethodId = "loadLabels";
   const char * const kLoadMappingMethodId = "loadMapping";
   const char * const kLoadTreeMethodId = "loadTree";
   const char * const kParseDataMethodId = "parseData";
   const char * const kComputeMethodId = "compute";
   const char * const kGetTreeMethodId = "getTree";
   const char * const kGetScoreMethodId = "getScore";
   const char * const kDestroyMethodId = "destroy";

   static const char kMessageArgumentSeparator = ':';

   template<class T> std::string toString(T x) {
     std::ostringstream str;
     str << x;
     return str.str();
   }

#ifdef INCLUDE_DATA
#include "data.h"
#endif
}  // namespace

// Note to the user: This glue code reflects the current state of affairs.  It
// may change.  In particular, interface elements marked as deprecated will
// disappear sometime in the near future and replaced with more elegant
// interfaces.  As of the time of this writing, the new interfaces are not
// available so we have to provide this code as it is written below.

/// The Instance class.  One of these exists for each instance of your NaCl
/// module on the web page.  The browser will ask the Module object to create
/// a new Instance for each occurence of the <embed> tag that has these
/// attributes:
///     type="application/x-nacl"
///     src="darwintree.nmf"
/// To communicate with the browser, you must override HandleMessage() for
/// receiving messages from the borwser, and use PostMessage() to send messages
/// back to the browser.  Note that this interface is entirely asynchronous.
class DarwintreeInstance : public pp::Instance {
   friend void* do_compute(void*);
   friend void* do_parse(void*);
   friend void do_delete_urlhandlers(void *user_data, int32_t result);

   public:
   /// The constructor creates the plugin-side instance.
   /// @param[in] instance the handle to the browser-side plugin instance.
   explicit DarwintreeInstance(PP_Instance instance) : pp::Instance(instance)
   {
      mstc = NULL;
      urlhandler[0] = urlhandler[1] = urlhandler[2] = urlhandler[3] = NULL;
      error = "";
      computation_running = false;
      parsing_running = false;
   }
   virtual ~DarwintreeInstance()
   {
      if(this->urlhandler[0]) delete this->urlhandler[0];
      if(this->urlhandler[1]) delete this->urlhandler[1];
      if(this->urlhandler[2]) delete this->urlhandler[2];
      if(this->urlhandler[3]) delete this->urlhandler[3];
      if(computation_running) pthread_kill(compute_thread,SIGKILL);
      if(parsing_running) pthread_kill(compute_thread,SIGKILL);
      if(this->mstc) delete this->mstc;
   }

   /// Handler for messages coming in from the browser via postMessage().  The
   /// @a var_message can contain anything: a JSON string; a string that encodes
   /// method names and arguments; etc.  For example, you could use
   /// JSON.stringify in the browser to create a message that contains a method
   /// name and some parameters, something like this:
   ///   var json_message = JSON.stringify({ "myMethod" : "3.14159" });
   ///   nacl_module.postMessage(json_message);
   /// On receipt of this message in @a var_message, you could parse the JSON to
   /// retrieve the method name, match it to a function call, and then call it
   /// with the parameter.
   /// @param[in] var_message The message posted by the browser.
   virtual void HandleMessage(const pp::Var& var_message) {
      if (!var_message.is_string()) {
         PostMessage(pp::Var("HandleMessage: FAILED (unknown message type)"));
      }
      std::string message = var_message.AsString();

      if(message == kParseDataMethodId) {
         this->parseData();
      } else if(message == kComputeMethodId) {
         this->compute();
      } else if(message == kGetTreeMethodId) {
         this->getTree();
      } else if(message == kGetScoreMethodId) {
         this->getScore();
      } else if(message == kDestroyMethodId) {
         this->destroy();
      } else if (message.find(kLoadMatricesMethodId) == 0) {
         size_t sep_pos = message.find_first_of(kMessageArgumentSeparator);
         if (sep_pos != std::string::npos) {
            std::string string_arg = message.substr(sep_pos + 1);
            this->loadMatrices(string_arg);
         } else {
            PostMessage(pp::Var("HandleMessage: FAILED (separator not found)"));
         }
      } else if (message.find(kLoadMappingMethodId) == 0) {
         size_t sep_pos = message.find_first_of(kMessageArgumentSeparator);
         if (sep_pos != std::string::npos) {
            std::string string_arg = message.substr(sep_pos + 1);
            this->loadMapping(string_arg);
         } else {
            PostMessage(pp::Var("HandleMessage: FAILED (separator not found)"));
         }
      } else if (message.find(kLoadLabelsMethodId) == 0) {
         size_t sep_pos = message.find_first_of(kMessageArgumentSeparator);
         if (sep_pos != std::string::npos) {
            std::string string_arg = message.substr(sep_pos + 1);
            this->loadLabels(string_arg);
         } else {
            PostMessage(pp::Var("HandleMessage: FAILED (separator not found)"));
         }
      } else if (message.find(kLoadTreeMethodId) == 0) {
         size_t sep_pos = message.find_first_of(kMessageArgumentSeparator);
         if (sep_pos != std::string::npos) {
            std::string string_arg = message.substr(sep_pos + 1);
            this->loadTree(string_arg);
         } else {
            PostMessage(pp::Var("HandleMessage: FAILED (separator not found)"));
         }
      } else {
         PostMessage(pp::Var("HandleMessage: FAILED (unknown method)"));
      }
   }

   void loadMatrices(const std::string &url);
   void loadMapping(const std::string &url);
   void loadLabels(const std::string &url);
   void loadTree(const std::string &url);
   void parseData();
   void compute();
   void getTree();
   void getScore();
   void destroy();

   void PostMessageOnMain(const pp::Var &msg);
   void DeleteUrlhandlers();

   private:
   UrlHandler* getUrl(const std::string &url, const std::string &name);

   MinSquareTreeCollection *mstc;
   std::string error;
   UrlHandler *urlhandler[4];
   pthread_t compute_thread;
   bool computation_running;
   bool parsing_running;
};

UrlHandler* DarwintreeInstance::getUrl(const std::string &url, const std::string &name) {
   UrlHandler *handler = new UrlHandler(this,name);
   std::cerr << "Requesting URL: " << url << std::endl;
   handler->processURL(url.c_str());
   return handler;
}

struct post_message_t {
   pp::Var message;
   DarwintreeInstance *instance;
};

static void do_post_message(void *user_data, int32_t result) {
   (void)result;
   post_message_t *pm = reinterpret_cast<post_message_t*>(user_data);
   pm->instance->PostMessage(pm->message);
   delete pm;
}

void DarwintreeInstance::PostMessageOnMain(const pp::Var &msg) {
   post_message_t *pm = new post_message_t;
   pm->instance = this;
   pm->message = msg;

   pp::Core *core = pp::Module::Get()->core();
   core->CallOnMainThread(0, pp::CompletionCallback(do_post_message, (void*)pm), 0);
}

void do_delete_urlhandlers(void *user_data, int32_t result) {
   (void)result;
   DarwintreeInstance *instance = reinterpret_cast<DarwintreeInstance*>(user_data);
   for(int i=0; i<4; ++i) {
      if(instance->urlhandler[i]) {
         delete instance->urlhandler[i];
         instance->urlhandler[i] = NULL;
      }
   }
}

void DarwintreeInstance::DeleteUrlhandlers() {
   pp::Core *core = pp::Module::Get()->core();
   core->CallOnMainThread(0, pp::CompletionCallback(do_delete_urlhandlers, (void*)this), 0);
}

void* do_parse(void* data) {
   std::vector<MinSquareTreeCollection::DblMatrix> matrices;
   MinSquareTreeCollection::IntMatrix mapping;
   std::vector<std::string> labels;
   PhyTree *tree = NULL;

   DarwintreeInstance *instance = reinterpret_cast<DarwintreeInstance*>(data);

#if !defined(INCLUDE_DATA)
   std::cerr << "started parsing" << std::endl;
   try {
      matrices = ProblemParser::parse_matrices(instance->urlhandler[0]->getData());
      //FIXME why does it segfault when I delete them here???
      //delete instance->urlhandler[0];
      //instance->urlhandler[0] = NULL;

      mapping = ProblemParser::parse_mapping(instance->urlhandler[1]->getData());
      //delete instance->urlhandler[1];
      //instance->urlhandler[1] = NULL;

      labels = ProblemParser::parse_labels(instance->urlhandler[2]->getData());
      //delete instance->urlhandler[2];
      //instance->urlhandler[2] = NULL;

      tree = ProblemParser::parse_tree(instance->urlhandler[3]->getData());
      //delete instance->urlhandler[3];
      //instance->urlhandler[3] = NULL;

      //FIXME ugly fix
      instance->DeleteUrlhandlers();

      std::cerr << "finished parsing" << std::endl;
   }
   catch(std::exception &e) {
      instance->error = e.what();
      instance->PostMessageOnMain(pp::Var(std::string(kParseDataMethodId) + ": FAILED (" + e.what() + ")"));

      //for(int i=0; i<4; ++i) {
      //   if(instance->urlhandler[i]) {
      //      delete instance->urlhandler[i];
      //      instance->urlhandler[i] = NULL;
      //   }
      //}

      if(tree) {
         delete tree;
      }
      instance->parsing_running = false;
      pthread_exit(NULL);
      return NULL;
   }

#else
   matrices = ProblemParser::parse_matrices(matrices_data);
   mapping = ProblemParser::parse_mapping(mapping_data);
   labels = ProblemParser::parse_labels(labels_data);
   tree = ProblemParser::parse_tree(tree_data);
#endif

   try {
      std::cerr << "initializing mstc" << std::endl;
      instance->mstc = new MinSquareTreeCollection(matrices,mapping,labels,*tree);
      instance->PostMessageOnMain(pp::Var(std::string(kParseDataMethodId) + ": OK"));
   }
   catch(std::exception &e) {
      if(instance->mstc) {
         delete instance->mstc;
         instance->mstc = NULL;
      }
      instance->PostMessageOnMain(pp::Var(std::string(kParseDataMethodId) + ": FAILED (" + e.what() + ")"));
   }
   if(tree) {
      delete tree;
   }

   instance->parsing_running = false;
   pthread_exit(NULL);
   return NULL;
}

void DarwintreeInstance::parseData() {
   if(parsing_running) {
      PostMessage(pp::Var(std::string(kParseDataMethodId) + ": FAILED (parsing already started)"));
      return;
   }

   if(computation_running) {
      PostMessage(pp::Var(std::string(kParseDataMethodId) + ": FAILED (computation running)"));
      return;
   }

   if(mstc) {
      delete mstc;
      mstc = NULL;
   }

#if !defined(INCLUDE_DATA)
   if(!urlhandler[0] || !urlhandler[0]->isFinished() || urlhandler[0]->hadError() ||
         !urlhandler[1] || !urlhandler[1]->isFinished() || urlhandler[1]->hadError() ||
         !urlhandler[2] || !urlhandler[2]->isFinished() || urlhandler[2]->hadError() ||
         !urlhandler[3] || !urlhandler[3]->isFinished() || urlhandler[3]->hadError()) {
      PostMessage(pp::Var(std::string(kParseDataMethodId) + ": FAILED (data not ready)"));
      return;
   }
#endif

   parsing_running = true;
   int ret = pthread_create(&compute_thread, NULL, &do_parse, (void*)this);
   if(ret != 0) {
      parsing_running = false;
      PostMessage(pp::Var(std::string(kParseDataMethodId) + ": FAILED (error starting compute thread)"));
   }
}

void DarwintreeInstance::loadMatrices(const std::string &url) {
   if(parsing_running) {
      PostMessage(pp::Var(std::string(kLoadMatricesMethodId) + ": FAILED (parsing running)"));
      return;
   }

#if !defined(INCLUDE_DATA)
  if(urlhandler[0]) {
    delete urlhandler[0];
    urlhandler[0] = NULL;
  }

   urlhandler[0] = getUrl(url,kLoadMatricesMethodId);
#else
   PostMessage(pp::Var(std::string(kLoadMatricesMethodId) + ": FAILED (data included)"));
#endif
}

void DarwintreeInstance::loadMapping(const std::string &url) {
   if(parsing_running) {
      PostMessage(pp::Var(std::string(kLoadMappingMethodId) + ": FAILED (parsing running)"));
      return;
   }

#if !defined(INCLUDE_DATA)
  if(urlhandler[1]) {
    delete urlhandler[1];
    urlhandler[1] = NULL;
  }

   urlhandler[1] = getUrl(url,kLoadMappingMethodId);
#else
   PostMessage(pp::Var(std::string(kLoadMappingMethodId) + ": FAILED (data included)"));
#endif
}

void DarwintreeInstance::loadLabels(const std::string &url) {
   if(parsing_running) {
      PostMessage(pp::Var(std::string(kLoadLabelsMethodId) + ": FAILED (parsing running)"));
      return;
   }

#if !defined(INCLUDE_DATA)
  if(urlhandler[2]) {
    delete urlhandler[2];
    urlhandler[2] = NULL;
  }

   urlhandler[2] = getUrl(url,kLoadLabelsMethodId);
#else
   PostMessage(pp::Var(std::string(kLoadLabelsMethodId) + ": FAILED (data included)"));
#endif
}

void DarwintreeInstance::loadTree(const std::string &url) {
   if(parsing_running) {
      PostMessage(pp::Var(std::string(kLoadTreeMethodId) + ": FAILED (parsing running)"));
      return;
   }

#if !defined(INCLUDE_DATA)
  if(urlhandler[3]) {
    delete urlhandler[3];
    urlhandler[3] = NULL;
  }

   urlhandler[3] = getUrl(url,kLoadTreeMethodId);
#else
   PostMessage(pp::Var(std::string(kLoadTreeMethodId) + ": FAILED (data included)"));
#endif
}

void* do_compute(void* data) {
   DarwintreeInstance *instance = reinterpret_cast<DarwintreeInstance*>(data);
   try {
      instance->mstc->compute(false);
      instance->PostMessageOnMain(pp::Var(std::string(kComputeMethodId) + ": OK"));
   }
   catch(std::exception &e) {
      instance->error = e.what();
      instance->PostMessageOnMain(pp::Var(std::string(kComputeMethodId) + ": FAILED (" + e.what() + ")"));
   }
   instance->computation_running = false;
   pthread_exit(NULL);
   return NULL;
}

void DarwintreeInstance::compute() {
  if(mstc) {
    if(computation_running || mstc->isComputed()) {
      PostMessage(pp::Var(std::string(kComputeMethodId) + ": FAILED (computation already started)"));
    } else if(mstc->isInitialized()) {
      computation_running = true;
      int ret = pthread_create(&compute_thread, NULL, &do_compute, (void*)this);
      if(ret != 0) {
        computation_running = false;
        PostMessage(pp::Var(std::string(kComputeMethodId) + ": FAILED (error starting compute thread)"));
      }
    } else {
      PostMessage(pp::Var(std::string(kComputeMethodId) + ": FAILED (not initialized yet)"));
    }
  } else {
    PostMessage(pp::Var(std::string(kComputeMethodId) + ": FAILED (not initialized yet)"));
  }
}

void DarwintreeInstance::getTree() {
  if(mstc) {
    if(mstc->isComputed()) {
      PostMessage(pp::Var(std::string(kGetTreeMethodId) + ": " + mstc->getTree()));
    } else {
      PostMessage(pp::Var(std::string(kGetTreeMethodId) + ": FAILED (result not computed yet)"));
    }
  } else {
    PostMessage(pp::Var(std::string(kGetTreeMethodId) + ": FAILED (result not computed yet)"));
  }
}

void DarwintreeInstance::getScore() {
  if(mstc) {
    if(mstc->isComputed()) {
      PostMessage(pp::Var(std::string(kGetScoreMethodId) + ": " + toString(mstc->getScore())));
    } else {
      PostMessage(pp::Var(std::string(kGetScoreMethodId) + ": FAILED (result not computed yet)"));
    }
  } else {
    PostMessage(pp::Var(std::string(kGetScoreMethodId) + ": FAILED (result not computed yet)"));
  }
}

void DarwintreeInstance::destroy() {
   if(computation_running) {
      pthread_kill(compute_thread,SIGKILL);
      computation_running = false;
   }
   if(parsing_running) {
      pthread_kill(compute_thread,SIGKILL);
      parsing_running = false;
   }
   if(mstc) {
      delete mstc;
      mstc = NULL;
   }
   for(int i=0; i<4; ++i) {
      if(urlhandler[i]) {
         delete urlhandler[i];
         urlhandler[i] = NULL;
      }
   }
}

/// The Module class.  The browser calls the CreateInstance() method to create
/// an instance of your NaCl module on the web page.  The browser creates a new
/// instance for each <embed> tag with type="application/x-nacl".
class DarwintreeModule : public pp::Module {
 public:
  DarwintreeModule() : pp::Module() {}
  virtual ~DarwintreeModule() {}

  /// Create and return a DarwintreeInstance object.
  /// @param[in] instance The browser-side instance.
  /// @return the plugin-side instance.
  virtual pp::Instance* CreateInstance(PP_Instance instance) {
    return new DarwintreeInstance(instance);
  }
};

namespace pp {
/// Factory function called by the browser when the module is first loaded.
/// The browser keeps a singleton of this module.  It calls the
/// CreateInstance() method on the object you return to make instances.  There
/// is one instance per <embed> tag on the page.  This is the main binding
/// point for your NaCl module with the browser.
Module* CreateModule() {
  return new DarwintreeModule();
}
}  // namespace pp
