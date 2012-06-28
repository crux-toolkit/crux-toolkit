/**
 * \file DelimitedFileWriter.cpp
 * $Revision: 1.0 $
 * DATE: October 19, 2010
 * AUTHOR: Barbara Frewen
 * \brief Object for writing tab-delimited files.
 * This class writes to files one line at a time, creating a table of
 * data.  The values for the current row are set one column at a time
 * and then the row is written to file, separating columns of data
 * with the delimiter.  The default delimiter is tab, but can be set
 * to any character.  A header row may be defined and printed to the
 * file at any row in the flile.  Once a header has been written,
 * every row after that will have the same number of fields.
 */

#include "DelimitedFileWriter.h"
#include "crux-file-utils.h"

using namespace std;

/**
 * \returns An empty DelimitedFileWriter object.
 */
DelimitedFileWriter::DelimitedFileWriter()
: file_ptr_(NULL),
  delimiter_('	') // default is tab
{
}

/**
 * \returns A DelimitedFileWriter object with the given file to
 * write to.
 */
DelimitedFileWriter::DelimitedFileWriter
(const char* filename) // full path of file
: file_ptr_(NULL),
  delimiter_('	') // default is tab
{
  this->openFile(filename);
}

/**
 * Destructor
 */
DelimitedFileWriter::~DelimitedFileWriter(){
  if( file_ptr_ ){
    file_ptr_->close();
    delete file_ptr_;
  }
}

/**
 * Writes any existing data, closes any open file and opens the
 * given file.
 */
void DelimitedFileWriter::openFile(const char* filename){
  // write any existing data and close file
  if( file_ptr_ ){
    writeRow();
    file_ptr_->close();
    delete file_ptr_;
  }

  // open the file if either it doesn't exist or if we are allowed to overwrite
  file_ptr_ = create_file(filename, get_boolean_parameter("overwrite"));
  if( file_ptr_ == NULL ){
    carp(CARP_FATAL, "Error creating file '%s'.", filename);
  }
}

/**
 * Sets the delimeter to separate columns.
 */
void DelimitedFileWriter::setDelimiter(char delimiter){
  delimiter_ = delimiter;
}

/**
 * Sets the name of the column at the given index, beginning with
 * zero.  Replaces any existing name.
 */
void DelimitedFileWriter::setColumnName(const string& name, ///< new name to set
                                        unsigned int col_idx){///< col to name
  while( column_names_.size() <= col_idx ){
    column_names_.push_back("");
  }
  column_names_[col_idx] = name;
}

/**
 * Sets the names of all columns, clearing any existing names.  Pass an
 * empty vector as the 'names' argument to clear any existing column
 * names.
 */
void DelimitedFileWriter::setColumnNames(const vector<string>& names){
  column_names_.clear();
  column_names_ = names;
}

/**
 * Writes the data of the current line to file, clears current data.
 * If there are column headers set, the line printed will have at
 * least as many fields as there are column headers.
 */
void DelimitedFileWriter::writeRow(){
  if( current_row_.empty() ){
    return;
  }

  // make the row as long as the header
  while( current_row_.size() < column_names_.size()){
    current_row_.push_back("");
  }
  // TODO? warning if row is longer than non-empty header?

  // print each value separated by delimiter
  *file_ptr_ << current_row_[0];
  for(size_t idx = 1; idx < current_row_.size(); idx++){
    *file_ptr_ << delimiter_ << current_row_[idx];
  }
  // end with newline
  *file_ptr_ << endl;

  // clear the current_row and refill with blanks
  // if there is a header, that is the min length
  // if not, each row can be a different length
  current_row_.assign(column_names_.size(), "");

}

/**
 * Writes any given column headers to file.  For missing column
 * headers (e.g. if columns 0 and 2 are set, but not 1) print
 * "column_#".
 */
void DelimitedFileWriter::writeHeader(){
  
  if( column_names_.empty() ){
    return;
  }
  
  if( file_ptr_ == NULL || !file_ptr_->is_open() ){
    carp(CARP_FATAL, "Cannot write to NULL delimited file.");
  }
  
  *file_ptr_ << column_names_[0];
  for(size_t idx = 1; idx < column_names_.size(); idx++){
    if( column_names_[idx].empty() ){
      *file_ptr_ << delimiter_ << "column_" << (idx+1); 
    } else {
      *file_ptr_ << delimiter_ << column_names_[idx]; 
    }
  }
  
  *file_ptr_ << endl;

  // with a header, each line must be that length
  current_row_.assign(column_names_.size(), "");
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

