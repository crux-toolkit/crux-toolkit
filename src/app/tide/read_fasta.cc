#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <gflags/gflags.h>
#include "records.h"
#include "raw_proteins.pb.h"
#include "abspath.h"

using namespace std;
using namespace google::protobuf::io;

#define CHECK(x) GOOGLE_CHECK(x)

class FASTAReader {
 public:
  FASTAReader(ZeroCopyInputStream* input) 
    : input_(input), id_(0) {
    eof_ = !ReadUntil('>',  &FASTAReader::DoNothing);
  }

  static void Init() {
    InitIsResidue();
  }

  const pb::Protein* GetProtein();

 private:
  typedef void (FASTAReader::*Action)(const char* buffer, int size); 
  bool ReadUntil(char c, Action what_to_do);

  void DoNothing(const char* buffer, int size) {}
  void AddToName(const char* buffer, int size) {
    *(current_protein_.mutable_name()) += string(buffer, size);
  }
  void AddToResidues(const char* buffer, int size) {
    vector<char> dest_buf(size);
    int total = 0;
    for (int i = 0; i < size; ++i, ++buffer)
      if (is_residue[*buffer])
        dest_buf[total++] = *buffer;
    *(current_protein_.mutable_residues()) += string(&(dest_buf[0]), total);
  }

  static bool is_residue[];
  static void InitIsResidue() {
    for (int i = 0; i < 256; ++i)
      is_residue[i] = false;
    const char* residues = "ACDEFGHIKLMNPQRSTVWY";
    for (const char* r = residues; *r; ++r)
      is_residue[*r] = true;
  }

  static string ShortenName(const string& name);

  ZeroCopyInputStream* input_;
  int id_;
  pb::Protein current_protein_;
  bool eof_;
};

bool FASTAReader::is_residue[256];

const pb::Protein* FASTAReader::GetProtein() {
  if (eof_)
    return NULL;
  current_protein_.Clear(); 
  // Assume we've read a '>'
  eof_ = !ReadUntil('\n', &FASTAReader::AddToName);
  if (!eof_)
    eof_ = !ReadUntil('>',  &FASTAReader::AddToResidues);
  current_protein_.set_id(id_++);
  current_protein_.set_name(ShortenName(current_protein_.name())); 
  return &current_protein_;
}


bool FASTAReader::ReadUntil(char c, Action what_to_do) {
  const char* buffer;
  int size;

  bool done = false;
  while (!done) {
    if (!input_->Next((const void**) &buffer, &size))
      return false;
    if  (size == 0) 
      continue;
    const char* pos = (const char*) memchr(buffer, c, size);
    if (pos != NULL) {
      int amt_used = pos - buffer;
      input_->BackUp(size - amt_used - 1); // -1 consumes separator
      size = amt_used;
      done = true;
    }
    (this->*what_to_do)(buffer, size);
  }
  return true;
}

string FASTAReader::ShortenName(const string& name) {
  const char* ws = " \t\v\n";
  int start = name.find_first_not_of(ws);
  if (start == string::npos)
    return "";
  int len = name.find_first_of(ws, start);
  return name.substr(start, len);
}

void TranslateFastaToPB(const string& fasta_filename,
                        const string& proteins_filename,
                        string* command_line = NULL,
			                  pb::Header* header = NULL) {
  int fd = open(fasta_filename.c_str(), O_RDONLY);
  CHECK(fd >= 0);
  FileInputStream input(fd);
  FASTAReader::Init();
  FASTAReader reader(&input);
  
  pb::Header tmp_header;
  if (header == NULL)
    header = &tmp_header;
  header->set_file_type(pb::Header::RAW_PROTEINS);
  if (NULL != command_line) 
    header->set_command_line(*command_line);
  pb::Header_Source* source = header->add_source();
  source->set_filename(AbsPath(fasta_filename));
  source->set_filetype("fasta");
  HeadedRecordWriter writer(proteins_filename, *header);
  CHECK(writer.OK());
  
  const pb::Protein* protein;
  while (protein = reader.GetProtein())
    writer.Write(protein);

  CHECK(writer.OK());
  close(fd);
}
