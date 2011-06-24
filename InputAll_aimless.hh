#ifndef INPUTALL_AIMLESS_HEADER
#define INPUTALL_AIMLESS_HEADER

#include "keywords_aimless.hh"
#include "Output.hh"


namespace phaser_io {
  class InputAll :
    public RESO,
    public PARTIALS,
    public EXCLUDE,
    public ANOMALOUS,
    public SCALES,
    public RUNSET,
    public BINS,
    public REJECT,
    public TIE,
    public NAME,
    public REFINE,
    public TITLE,
    public ONLYMERGE,
    public BLANK,
    public SDCORRECTION,
    public INTENSITIES,
    public KEEP,
    public OUTPUT,
    public DUMP,
    public RESTORE,
    public ANALYSIS

  {
  public:
    InputAll(Preprocessor&); 
    InputAll();
    // Entry to do nothing if "online"  ie if stdin not connected to file)
    InputAll(const bool& OnLine);
    InputAll(const bool& OnLine, Output& output);  // ... with echo to output
    ~InputAll();
    void Analyse();

    // true is input has been read (ie not online)
    bool Valid() const {return !(OnLine_);}

  private:
    bool OnLine_;

    void init(const bool& OnLine, const bool& echo, Output& output);

  };
    
} // phaser_io
  
#endif

