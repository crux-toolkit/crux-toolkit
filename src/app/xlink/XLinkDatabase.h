#ifndef XLINKDATABASE_H_
#define XLINKDATABASE_H_
#include "model/objects.h"
#include "model/Peptide.h"
#include "model/Database.h"
#include "XLinkablePeptide.h"
#include "SelfLoopPeptide.h"
#include "LinearPeptide.h"
#include "MonoLinkPeptide.h"

#include <vector>

class XLinkDatabase {

 protected:
  static Database* protein_database_;
  static XLinkBondMap bondmap_;


  static std::vector<std::vector<Crux::Peptide*> > target_peptides_; 
    ///< all target peptides generated, indexed by missed cleavages

  static std::vector<std::vector<Crux::Peptide*> > decoy_peptides_; 
    ///< all decoy peptides generated, indexed by missed cleavages

  static std::vector<LinearPeptide> target_linear_peptides_;
  static std::vector<LinearPeptide> decoy_linear_peptides_;

  static std::vector<MonoLinkPeptide> target_monolink_peptides_;

  static std::vector<SelfLoopPeptide> target_selfloop_peptides_;
  static std::vector<SelfLoopPeptide> decoy_selfloop_peptides_;

  static std::vector<XLinkablePeptide> target_xlinkable_peptides_;
  static std::vector<XLinkablePeptide> decoy_xlinkable_peptides_;

  static std::vector<XLinkablePeptide> target_xlinkable_peptides_flatten_;
  static std::vector<XLinkablePeptide> target_xlinkable_peptides2_; //Peptides that could be selfloops.
  static std::vector<XLinkablePeptide> decoy_xlinkable_peptides2_;  

  static std::vector<XLinkablePeptide> decoy_xlinkable_peptides_flatten_;

  static void findLinearPeptides(
    vector<Crux::Peptide*>& peptides, 
    vector<LinearPeptide>& linears
  );
  static void generateAllLinears(bool decoy);
  static void generateAllLinkablePeptides(
    std::vector<Crux::Peptide*>& peptides, 
    std::vector<XLinkablePeptide>& xpeptides);
  
  static void findSelfLoops(
    std::vector<XLinkablePeptide>& linkable_peptides,
    std::vector<SelfLoopPeptide>& ans);

  static void generateAllSelfLoops(bool decoy);
  static void flattenLinkablePeptides(
   std::vector<XLinkablePeptide>& xpeptides,
   std::vector<XLinkablePeptide>& flattened
   );

  static void filterLinkablePeptides(
    std::vector<XLinkablePeptide>& xpeptides,
    std::vector<XLinkablePeptide>& filtered_xpeptides
    );  

  static bool addPeptideToDatabase(Crux::Peptide* peptide);  
    
 public:
  XLinkDatabase() {;}
  virtual ~XLinkDatabase() {;}

  static void initialize();
  static void finalize();
  
  static XLinkBondMap& getXLinkBondMap();

  static std::vector<XLinkablePeptide>::iterator getXLinkableBegin();
  static std::vector<XLinkablePeptide>::iterator getXLinkableBegin(FLOAT_T min_mass);
  static std::vector<XLinkablePeptide>::iterator getXLinkableEnd();
  static int getNLinkable();
  static std::vector<XLinkablePeptide>& getXLinkablePeptides(
    bool decoy
  );

  static XLinkablePeptide& getXLinkablePeptide(
    bool decoy, 
    int idx
  );

  static std::vector<XLinkablePeptide>::iterator getXLinkableFlattenBegin();
  static std::vector<XLinkablePeptide>::iterator getXLinkableFlattenBegin(
    bool decoy, 
    FLOAT_T min_mass
  );
  static std::vector<XLinkablePeptide>::iterator getXLinkableFlattenEnd();
  
  static std::vector<XLinkablePeptide>::iterator getXLinkableFlattenEnd(
    bool decoy,
    FLOAT_T max_mass
  );

  static std::vector<SelfLoopPeptide>::iterator getSelfLoopBegin(
    bool decoy
  );

  static std::vector<SelfLoopPeptide>::iterator getSelfLoopEnd(
    bool decoy
  );
  
  static std::vector<SelfLoopPeptide>::iterator getSelfLoopBegin(
    bool decoy, 
    FLOAT_T min_mass
  );

  static std::vector<LinearPeptide>::iterator getLinearBegin(
    bool decoy
  );
  
  static std::vector<LinearPeptide>::iterator getLinearBegin(
    bool decoy,
    FLOAT_T min_mass
  );

  static std::vector<LinearPeptide>::iterator getLinearEnd(
    bool decoy
  );

  static std::vector<LinearPeptide>::iterator getLinearEnd(
    bool decoy,
    FLOAT_T max_mass
  );
  
  static std::vector<LinearPeptide>::iterator getLinearEnd(
    bool decoy,
    std::vector<LinearPeptide>::iterator &siter,
    FLOAT_T max_mass
  );


  static std::vector<MonoLinkPeptide>::iterator getMonoLinkBegin(
    bool decoy
  );
  
  static std::vector<MonoLinkPeptide>::iterator getMonoLinkBegin(
    bool decoy,
    FLOAT_T min_mass
  );

  static std::vector<MonoLinkPeptide>::iterator getMonoLinkEnd(
    bool decoy
  );

  static std::vector<MonoLinkPeptide>::iterator getMonoLinkEnd(
    bool decoy,
    FLOAT_T max_mass
  );
  
  static std::vector<MonoLinkPeptide>::iterator getMonoLinkEnd(
    bool decoy,
    std::vector<MonoLinkPeptide>::iterator &siter,
    FLOAT_T max_mass
  );


  
  //static std::vector<std::pair<int, vector<XLinkablePeptide> > >& getTargetProteinIdxToXPeptides();


  static void print();

};
#endif
