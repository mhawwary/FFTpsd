

#include"general_tools.h"

// Sorting double array
void QuickSort(double*& szArray,double*& ssArray , int nLower, int nUpper)
{
  // Check for non-base case
  if (nLower < nUpper)
  {
    // Split and sort partitions
    int nSplit = Partition(szArray,ssArray, nLower, nUpper);
    QuickSort(szArray,ssArray, nLower, nSplit - 1);
    QuickSort(szArray,ssArray, nSplit + 1, nUpper);
  }
}

// QuickSort partition implementation
int Partition(double*& szArray,double*& ssArray, int nLower, int nUpper)
{
  // Pivot with first element
  int nLeft = nLower + 1;
  double szPivot = szArray[nLower];
  //double ssPivot=ssArray[nLower];
  int nRight = nUpper;

  // Partition array elements
  double szSwap,ssSwap;
  while (nLeft <= nRight)
  {
    // Find item out of place
    while (nLeft <= nRight && szArray[nLeft]<= szPivot)
      nLeft = nLeft + 1;
    while (nLeft <= nRight && szArray[nRight]>szPivot)
      nRight = nRight - 1;

    // Swap values if necessary
    if (nLeft < nRight)
    {
      szSwap = szArray[nLeft];             ssSwap=ssArray[nLeft];
      szArray[nLeft] = szArray[nRight];    ssArray[nLeft]=ssArray[nRight];
      szArray[nRight] = szSwap;            ssArray[nRight] = ssSwap;
      nLeft = nLeft + 1;
      nRight = nRight - 1;
    }
  }

  // Move pivot element
  szSwap = szArray[nLower];
  szArray[nLower] = szArray[nRight];
  szArray[nRight] = szSwap;

  ssSwap = ssArray[nLower];
  ssArray[nLower] = ssArray[nRight];
  ssArray[nRight] = ssSwap;

  return nRight;
}

// Sorting two double arrays and an ID int one
void QuickSort3(double*& mainArray,double*& A1 , int*& A2 , int nLower, int nUpper)
{
  // Check for non-base case
  if (nLower < nUpper)
  {
    // Split and sort partitions
    int nSplit = Partition3(mainArray,A1,A2,nLower, nUpper);
    QuickSort3(mainArray,A1,A2, nLower, nSplit - 1);
    QuickSort3(mainArray,A1,A2, nSplit + 1, nUpper);
  }
}

// QuickSort partition implementation for three arrays
int Partition3(double*& mainArray,double*& A1 , int*& A2, int nLower, int nUpper)
{
  // Pivot with first element
  int nLeft = nLower + 1;
  double mainPivot =  mainArray[nLower];
  int nRight = nUpper;

  // Partition array elements
  double mainSwap,A1Swap;  int A2Swap;
  while (nLeft <= nRight)
  {
    // Find item out of place
    while (nLeft <= nRight && mainArray[nLeft]<= mainPivot)
      nLeft = nLeft + 1;
    while (nLeft <= nRight && mainArray[nRight]>mainPivot)
      nRight = nRight - 1;

    // Swap values if necessary
    if (nLeft < nRight)
    {
      mainSwap = mainArray[nLeft];     mainArray[nLeft] = mainArray[nRight];     mainArray[nRight] = mainSwap;

      A1Swap= A1[nLeft];  A1[nLeft]=A1[nRight];   A1[nRight] = A1Swap;
      A2Swap= A2[nLeft];  A2[nLeft]=A2[nRight];   A2[nRight] = A2Swap;

      nLeft = nLeft + 1;
      nRight = nRight - 1;
    }
  }

  // Move pivot element
  mainSwap = mainArray[nLower];   mainArray[nLower] = mainArray[nRight];    mainArray[nRight] = mainSwap;

  A1Swap = A1[nLower];    A1[nLower] = A1[nRight];    A1[nRight] = A1Swap;
  A2Swap = A2[nLower];    A2[nLower] = A2[nRight];    A2[nRight] = A2Swap;

  return nRight;
}

void open_inputfile_forreading(std::string& o_readfname_, std::ifstream& o_input_){
  o_input_.open(o_readfname_);
  if(!o_input_.is_open()){
    std::string o_error_message_;
    o_error_message_ = "Failed to open "+ o_readfname_ +", or file does not exist!";
    FatalError(o_error_message_.c_str());
    exit(0);
  }
}

std::string remove_extension(const std::string& filename){
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos) return filename;
  return filename.substr(0, lastdot);
}

std::string remove_from_end_up_to_string(const std::string field,const std::string& filename){
  size_t lastdot = filename.find_last_of(field);
  if (lastdot == std::string::npos) return filename;
  return filename.substr(0, lastdot);
}

std::string GetFileExtension(const std::string& FileName)
{
  if(FileName.find_last_of(".") != std::string::npos)
    return FileName.substr(FileName.find_last_of(".")+1);
  return "";
}


std::string GetFileDirectory(const std::string& FileName){
  size_t lastdot = FileName.find_last_of("/");
  if (lastdot == std::string::npos) return "./";
  std::string directory=FileName.substr(0, lastdot)+"/";
  return directory;
}

bool is_a_comment_line(const std::string& line_in){
  bool check=false;
  if(line_in.find_first_of("#") !=std::string::npos  ||
     line_in.find_first_of("%") !=std::string::npos  ||
     line_in.find_first_of("//")!=std::string::npos  ||
     line_in.find_first_of("/!")!=std::string::npos  ||
     line_in.find_first_of("/*")!=std::string::npos )
    check=true;
  return check;
}

bool is_a_text_line(const std::string& line_in){
  bool check=false;
  if(line_in.find_first_of("auoi") !=std::string::npos)
    check=true;
  else if(line_in.find_first_of("e") !=std::string::npos){
    std::size_t pos0=line_in.find_first_of("e");
    std::string str0=line_in.substr(pos0+1,pos0+3);
    if(is_a_text_line(str0))
      check=true;
    else if(is_a_text_line(line_in.substr(pos0-2,pos0)))
      check=true;
  }

  return check;
}

// This function converts line string of numbers to integers
void line2intdata (const std::string line, std::vector<int> &data, int& n_data_oneline)
{
  std::istringstream istr(line);
  int j = 0;
  std::string temp;
  char* c=nullptr;
  c = new char[line.size()+1];
  int tempint;

  while(!istr.eof()){
    istr >> temp;
    std::strcpy(c, temp.c_str());
    tempint = atoi(c);
    data[j] = tempint;
    j++;
  }

  n_data_oneline = j;

  emptyarray(c);

  return ;
}

// This function converts line string of numbers to integers
void line2doubledata (const std::string line, std::vector<double> &data, int& n_data_oneline)
{
  std::istringstream istr(line);
  int j = 0;
  double value;

  do{
    istr >> value;  // need to put something to skip the line if it contains letters or #
    data.push_back(value);
    j++;
  }while(!istr.eof() && j<1e5);

  n_data_oneline = j;

  return ;
}

