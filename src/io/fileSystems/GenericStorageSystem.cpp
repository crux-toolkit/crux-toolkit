#include "io/fileSystems/GenericStorageSystem.h"
#include "io/fileSystems/BoostFileSystem.h"
#include <boost/iostreams/filtering_stream.hpp>

#ifdef AWS      // build in AWS support only if the library is present at compilation time.
#include "io/fileSystems/AwsS3System.h"
#endif

#include <boost/algorithm/string.hpp>
#include "io/carp.h"
#include <string>
#include <regex>
#include <vector>
#include <thread>

//TODO_RC: use a map instead of plain array
GenericStorageSystem * GenericStorageSystem::m_system[] = {nullptr, nullptr};
std::map<int, GenericStorageSystem::StreamRecord> GenericStorageSystem::m_openStreams; 
std::mutex GenericStorageSystem::m_registerMutex;
int GenericStorageSystem::Identifiable::globalStreamCounter{0};

GenericStorageSystem* GenericStorageSystem::getStorage(string p_path)
{
    int storageIndex = getStorageIndex(p_path);
    if(m_system[storageIndex] == nullptr)  
        switch(storageIndex){
            case 0:
                m_system[storageIndex] = new BoostFileSystem();
                carp(CARP_INFO, "Using file:// storage system for path %s.", p_path.c_str());
                break;
#ifdef AWS                // build in AWS support only if the library is present at compilation time.
            case 1:
                m_system[storageIndex] = new AwsS3System();
                carp(CARP_INFO, "Using s3:// storage system for path %s.", p_path.c_str());
                break;
#endif                
            default:
                carp(CARP_FATAL, "Path %s is not supported", p_path.c_str());
        }
    return m_system[storageIndex];
}

void GenericStorageSystem::cleanup(){
    for(int i; i < MAX_STORAGES; i++){
        if(m_system[i] != nullptr)
            delete m_system[i];
    }
}

/*
* Private method. Determines the protocol specified in the path string and returns the 
* index of appropriate storage system.
*/
int GenericStorageSystem::getStorageIndex(const string &path){
    regex pattern {"^([a-z0-9]+):"};
    string protocol_name{"file"};
    vector<string> known_protocols{"file", "s3"};
    string lower_path = boost::algorithm::to_lower_copy(path);
    smatch matches;
    //extract protocal name from the URI path
    if(regex_search(lower_path, matches, pattern)){
        protocol_name = matches[1];
    }
    //find the protocol position
    int pos = 0;
    for(auto p = known_protocols.begin(); p != known_protocols.end(); ++p){
        if(protocol_name == *p)
            break;
        else{
            pos++;
        }
    }
    if(pos >= MAX_STORAGES)
        carp(CARP_FATAL, "Unknown storage system %s in resource URI %s", protocol_name.c_str(), path.c_str());

    return pos;
}

void GenericStorageSystem::CloseStream(ios_base& stream){

        Identifiable* ident = dynamic_cast<Identifiable*>(&stream);
        if(ident == nullptr)
            carp(CARP_FATAL, "Cannot close a stream that does not have stream ID.");

        std::lock_guard<std::mutex> lg{m_registerMutex};

        auto pRec = m_openStreams.find(ident->streamId);
        if(pRec == m_openStreams.end()){
            carp(CARP_WARNING, "The stream with ID %d is not yet open", ident->streamId);
            return;
        carp(CARP_INFO, "Closing %s.", pRec->second.toString().c_str());
        m_system[pRec->second.system_id]->_CloseStream(pRec->second);
        m_openStreams.erase(ident->streamId);  
    }
}

void GenericStorageSystem::_RegisterStream(const GenericStorageSystem::StreamRecord& p_streamRec){
    if(m_openStreams.find(p_streamRec.stream_id) != m_openStreams.end()) 
    {
            carp(CARP_WARNING, "The stream with ID is already open.", p_streamRec.stream_id); 
            return ;
    }

    carp(CARP_INFO, "Registering %d", p_streamRec.toString().c_str());
    {
        std::lock_guard<std::mutex> lg{m_registerMutex};
        m_openStreams.insert({p_streamRec.stream_id, p_streamRec});
    }
}

istream& GenericStorageSystem::GetReadStream(const string &path){
    StreamRecord rec{m_systemId, true, 0, (void*)nullptr, path};
    istream* str = getReadStreamImpl(path, rec);
    IdentifiableInputStream* result = new IdentifiableInputStream(*str);
    rec.stream_id = result->streamId;
    _RegisterStream(rec);
    return *result;
  }



//if paths are different it opens input and output stream and 
// moves data over. Otherwise it delegates to CopyLocal
void GenericStorageSystem::Copy(const string &orig, const string &dest){

}
