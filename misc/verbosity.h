#ifndef __VERBOSITY_H__
#define __VERBOSITY_H__

// log4cpp has following levels of logging :
// 7 - debug
// 6 - info
// 5 - notice
// 4 - warn
// 3 - error
// 2 - critical
// 1 - alert
// 0 - fatal
// Verbosity v means that all messages <=v can be printed
// e.g. for verbosity 7, all messages shd be printed
// e.g. for verbosity 3, only error, critical, alert and fatal messages shd be printed

class Verbosity {
public :
    static Verbosity & instance();
    bool isLevelLEq(int level);
    void setVerbo(int level);
private :
    Verbosity();
    Verbosity(const Verbosity &);
    int verbo;
};

void setVerbosity(int v);

bool verbose(int level); // is verbosity at least at this level? ie is level <= verbosity

#endif // __VERBOSITY_H__
