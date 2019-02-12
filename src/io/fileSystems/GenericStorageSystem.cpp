#include "io/fileSystems/GenericStorageSystem.h"
#include "io/fileSystems/BoostFileSystem.h"
#include "io/fileSystems/AwsS3System.h"
#include <boost/algorithm/string.hpp>
#include "io/carp.h"
#include <string>
#include <regex>

//TODO_RC: use a map instead of plain array
GenericStorageSystem * GenericStorageSystem::m_system[] = {nullptr, nullptr};

GenericStorageSystem* GenericStorageSystem::getStorage(string p_path)
{
    int storageIndex = getStorageIndex(p_path);
    if(m_system[storageIndex] == nullptr)  
        switch(storageIndex){
            case 0:
                m_system[storageIndex] = new BoostFileSystem();
                carp(CARP_INFO, "Using file:// storage system for path %s.", p_path);
                break;
#ifdef AWS                // build in AWS support only if the library is present at compilation time.
            case 1:
                m_system[storageIndex] = new AwsS3System();
                carp(CARP_INFO, "Using s3:// storage system for path %s.", p_path);
                break;
#endif                
            default:
                carp(CARP_FATAL, "Path %s is not supported", p_path);
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
        carp(CARP_FATAL, "Unknown storage system %s in resource URI %s", protocol_name, path);

    return pos;
}

//if paths are different it opens input and output stream and 
// moves data over. Otherwise it delegates to CopyLocal
void GenericStorageSystem::Copy(const string &orig, const string &dest){

}