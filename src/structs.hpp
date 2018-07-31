/**
 * Structs used to hold information from GTF and SAM files.
 */

#ifndef structs_hpp
#define structs_hpp

#include <string>
#include <fstream>
#include <seqan/gff_io.h>
#include "util.hpp"

#define PARSE_FAILED 0xFFFF  // Read failed if read.pos is this number.

/**
 * Corresponds to an exon.
 */
struct Exon {
    int start;                        /* Start index of exon.
                                       *  1-indexed. */
    int end;                          /* End index of exon. 1-indexed. */
    std::string seqname;              /* Chromosome/scaffold of exon. */
    std::vector<int> *transcripts;    /* All transcripts that contain this
                                       * exon. Likely to contain multiple
                                       * entries, since exons are
                                       *  repetitive. */
    
    /**
     * Constructor for an exon from init_start to init_end.
     * @param init_start        start
     * @param init_end          end
     */
    Exon(uint init_start = 0, uint init_end = 0) {
        start = init_start;
        end = init_end;
        transcripts = new std::vector<int>;
    }
    
    /**
     * Constructor for an exon using info in a sequence seq. Will fill in all
     * fields with relevant information, not including transcripts.
     * @param seq       Sequence from which to copy information.
     */
    Exon(seqan::GffRecord &rec) {
        start = rec.beginPos;
        end = rec.endPos;
        seqname = lower(seqan::toCString(rec.ref));
        transcripts = new std::vector<int>;
    }

    /**
     * Copy constructor.
     * @param e         Exon to copy. const because if it's not, I get a
     * strange error. If Microsoft documentation can be this sketchy, mine can
     * be too.
     */
    Exon(const Exon &e) {
        start = e.start;
        end = e.end;
        seqname = e.seqname;
        transcripts = new std::vector<int>(*e.transcripts);
    }
    
    /**
     * Destructor.
     */
    ~Exon() {
        delete transcripts;
    }
};

#endif

