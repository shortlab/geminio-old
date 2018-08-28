#include "GeminioUtils.h"

#include <cstdlib>
#include <cctype>

namespace GeminioUtils
{

int
getGroupNumber(std::string str)
{
  int len = str.length(), i = len;

  while (std::isdigit(str[i-1]))
    i--;

  int no = std::atoi((str.substr(i)).c_str());
  while (i-- >= 0)
  {
    if (str[i] == 'v')
      return no;
    if (str[i] == 'i')
      return -no;
  }

  return no;
}

}
