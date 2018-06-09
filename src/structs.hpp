/**
 * Structs used to hold information from GTF and SAM files.
 */

#include <string>
#include <fstream>

#ifndef structs_hpp
#define structs_hpp

#define PARSE_FAILED 0xFFFF  // Read failed if read.pos is this number.

/**
 * Corresponds to a sequence in a GTF file. See GTF documentation for more
 * detail on each member.
 */
struct Sequence {
    uint start;          /*< Start index of sequence. 1-indexed. */
    uint end;            /*< End index of sequence. 1-indexed. */
    std::string seqname; /*< Chromosome or scaffold of sequence. */
    std::string feature; /*< Feature type, e.g. exon, transcript. */
    std::string id;      /*< Transcript ID, as indicated in attribute field. */
};

/**
 * Corresponds to an exon.
 */
struct Exon {
    uint32_t start;                        /*< Start index of exon.
                                            1-indexed. */
    uint32_t end;                          /*< End index of exon. 1-indexed. */
    std::string seqname;                   /*< Chromosome/scaffold of exon. */
    std::vector<uint64_t> *transcripts;    /*< All transcripts that contain this
                                            exon. Likely to contain multiple
                                            entries, since exons are
                                            repetitive. */
    
    /**
     * Constructor for an exon from init_start to init_end.
     * @param init_start        start
     * @param init_end          end
     */
    Exon(uint init_start = 0, uint init_end = 0) {
        start = init_start;
        end = init_end;
        transcripts = new std::vector<uint64_t>;
    }
    
    /**
     * Constructor for an exon using info in a sequence seq. Will fill in all
     * fields with relevant information, not including transcripts.
     * @param seq       Sequence from which to copy information.
     */
    Exon(Sequence &seq) {
        start = seq.start;
        end = seq.end;
        seqname = seq.seqname;
        transcripts = new std::vector<uint64_t>;
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
        transcripts = new std::vector<uint64_t>(*e.transcripts);
    }
    
    /**
     * Destructor.
     */
    ~Exon() {
        delete transcripts;
    }
};

/**
 * Corresponds to an entry for a read in a SAM file. This is not necessarily
 * one entire read, since a single (multi-mapping) read may have multiple
 * entries. See SAM documentation for more detail on each member.
 */
struct Read {
    int flag;           /*< Flag indicating some special status of entry. */
    uint32_t pos;           /*< Start index of entry. 1-indexed. */
    uint32_t end;           /*< End index of entry. Not necessarily equal to
                         pos + (read length), since base pairs may have been
                         deleted, spliced, etc. */
    std::string qname;  /*< Chromosome or scaffold of entry. */
    std::string rname;  /*< Name of read represented by this entry. */
    std::string cigar;  /*< CIGAR string. See documentation for details. */
};

#endif
