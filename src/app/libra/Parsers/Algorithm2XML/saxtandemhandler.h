#if !defined(TANSAXHANDLER_H)
#define TANSAXHANDLER_H

#include <expat.h>
#include <string>
#include <string.h>

using namespace std;

/**
* eXpat SAX parser wrapper.
*/
class TANSAXHandler
{
public:
	TANSAXHandler();
	virtual ~TANSAXHandler();

	/**
	* Receive notification of the start of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the start of
	* each element (such as allocating a new tree node or writing
	* output to a file).</p>
	*/
	virtual void startElement(const XML_Char *el, const XML_Char **attr);

	/**
	* Receive notification of the end of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the end of
	* each element (such as finalising a tree node or writing
	* output to a file).</p>
	*/
	virtual void endElement(const XML_Char *el);

	/**
	* Receive notification of character data inside an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method to take specific actions for each chunk of character data
	* (such as adding the data to a node or buffer, or printing it to
	* a file).</p>
	*/
	virtual void characters(const XML_Char *s, int len);

	/**
	* Open file and stream data to the SAX parser.  Must call
	* setFileName before calling this function.
	*/
	bool parse();

	inline void setFileName(const char* fileName)
	{
		m_strFileName = fileName;
	}

	// Helper functions
	inline bool isElement(const char *n1, const XML_Char *n2)
	{	return (strcmp(n1, n2) == 0); }

	inline bool isAttr(const char *n1, const XML_Char *n2)
	{	return (strcmp(n1, n2) == 0); }

	inline const char* getAttrValue(const char* name, const XML_Char **attr)
	{
		for (int i = 0; attr[i]; i += 2)
		{
			if (isAttr(name, attr[i]))
				return attr[i + 1];
		}

		return "";
	}

protected:
	XML_Parser m_parser;

	string  m_strFileName;
};

#endif              //TANSAXHANDLER_H
