#ifndef __SSTWEKAERED_H__
#define __SSTWEKAERED_H__

class SStweakerED {
public:
    SStweakerED(vector<int>& opInds, const char* desc);
    SStweakerED() {};
    virtual ~SStweakerED() {};

    void build(VVF & pts);
    vector<int>& getOP();
    float energy();
    const char* name();
protected:
    vector<int> op;
    string description;
};

#endif //__SSTWEKAERED_H__
