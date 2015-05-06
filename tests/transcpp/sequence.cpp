#include "sequence.h"


Sequence::Sequence() {}

Sequence::Sequence(const string& s)      { sequence = string2int(s); }

Sequence::Sequence(const vector<int>& s) { sequence = s; }

Sequence::Sequence(const vector<char>& s) { sequence = char2int(s); }



void Sequence::setSequence(const string& s)      { sequence = string2int(s); }
void Sequence::setSequence(const vector<int>& s) { sequence = s; }


vector<int>&  Sequence::getSequence() {return sequence; }
string Sequence::getSequenceString() {return int2string(sequence); }
