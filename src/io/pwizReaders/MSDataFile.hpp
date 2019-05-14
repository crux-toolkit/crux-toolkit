#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/Reader.hpp"
#include "io/fileSystems/GenericStorageSystem.h"
#include "io/pwizReaders/DefaultReaderList.hpp"

#include "io/carp.h"


using namespace std;

namespace crux {
namespace msdata {

    shared_ptr<DefaultReaderList> defaultReaderList_;

    struct MsDataFileAws : public pwiz::msdata::MSDataFile {

        //overwriting the parent class constructor to implement AWS reads
        MsDataFileAws(const std::string& filename, 
                const pwiz::msdata::Reader* reader = 0,
                bool calculateSourceFileChecksum = false)
        {
            // peek at head of file 
            //string head = read_file_header(filename, 512);
            GenericStorageSystem* s = GenericStorageSystem::getStorage(filename);
            if(!s-> IsRegularFile(filename)){
                carp(CARP_FATAL, "Path %s does not exist.", filename.c_str());
            }

            string head = s->Read(filename, 512);

            const pwiz::msdata::Reader* useThisReader;

            if (reader)
                useThisReader = reader;
            else{
                if (!defaultReaderList_.get())
                    defaultReaderList_ = shared_ptr<DefaultReaderList>(new DefaultReaderList());
                useThisReader = defaultReaderList_.get();
            }

            if (!useThisReader->accept(filename, head))
                carp(CARP_FATAL, "File %s format is invalid for reader %s.", filename.c_str(), useThisReader->getType());

            useThisReader->read(filename, head, *this);


            if (calculateSourceFileChecksum && !fileDescription.sourceFilePtrs.empty())
                calculateSourceFileSHA1(*fileDescription.sourceFilePtrs.back());
        }
    };
}
}