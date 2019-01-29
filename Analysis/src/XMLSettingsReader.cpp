/** @file XMLSettingsReader.cpp
 *  @brief Implementation of XMLSettingsReader.
 *
 *  Function definitions for XMLSettingsReader are provided.
 *  This class is the main  class to read xml configuration files for ZDCs and RPD.
 *  DOM logic is chosen for the xml implementation
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */


#include "XMLSettingsReader.h"

using namespace xercesc;

/** @brief Default Constructor for XMLSettingsReader.
 */
XMLSettingsReader::XMLSettingsReader( void ){

    try {
        XMLPlatformUtils::Initialize();
      } catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());

        std::cout << "Fatal error during the parser initialization! " << std::endl;
        XMLString::release(&message);
      }

      //Initial number of root nodes = 0
      rootNode = 0;

      //Initializing the DOM parser and the error handles
      myDOMParser = new XercesDOMParser();
      myDOMParser->setValidationScheme(XercesDOMParser::Val_Always);
      myDOMParser->setDoNamespaces(false);
      ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
      myDOMParser->setErrorHandler(errHandler);
      myDOMParser->setDoSchema(false);
      myDOMParser->setValidationConstraintFatal(false);

}

/** @brief Destructor for XMLSettingsReader.
 */
XMLSettingsReader::~XMLSettingsReader( void ){

    if (myDOMParser != NULL) delete myDOMParser;
    XMLPlatformUtils::Terminate();

}

/*! @brief This function returns the number of baseNodes with the specified name
 *  @param _nodeName the name of the baseNodes to count
 *
 * Assume:
 * <root>
 *  <detector>
 *   <channel></channel>
 *  </detector>
 *  <detector>
 *   <channel></channel>
 *   <channel></channel>
 *  </detector>
 * </root>
 * getBaseNodeCount("detector") returns 2
 *
 */

int XMLSettingsReader::getBaseNodeCount(std::string _nodeName){

    return rootNode->getElementsByTagName(XMLString::transcode(_nodeName.c_str()))->getLength();

}


/*! @brief returns the number of child Nodes in a give parent Node
 *  @param _nodeName parent node name
 *  @param _nodeNumber number of the parentNode
 *  @param _childName the name of the childNode to count
 *
 * Assume:
 * <root>
 *  <detector>
 *   <channel></channel>
 *  </detector>
 *  <detector>
 *   <channel></channel>
 *   <channel></channel>
 *  </detector>
 * </root>
 * getNodeChildCount("detector",0,"channel") returns 1
 * getNodeChildCount("detector",1,"channel") returns 2
 */

int XMLSettingsReader::getNodeChildCount(std::string _nodeName, int _nodeNumber, std::string _childName){

      XMLCh* bufferNode = XMLString::transcode(_nodeName.c_str());
      DOMNodeList* list = myDOMParser->getDocument()->getElementsByTagName(bufferNode);
      XMLString::release(&bufferNode);

      if (list->getLength() <= (unsigned int) _nodeNumber)
        return 0;

      DOMElement* parent = dynamic_cast<DOMElement*>(list->item(_nodeNumber));
      //Now we create a child list
      DOMNodeList* childList = parent->getElementsByTagName(XMLString::transcode(_childName.c_str()));
      //And we return its size
      return (int) childList->getLength();

}

/*! @brief This method returns a string the value of a child
 *  @param _baseNode the name of the baseNode
 *  @param _baseNumber the number of the baseNode to choose
 *  @param _childNode the name of the childNode to choose
 *
 * Usage example:
 * Assume:
 * <root>
 *  <detector>
 *   <channel>C1</channel>
 *   <delay>100</delay>
 *  </detector>
 *  <detector>
 *   <channel>C2</channel>
 *   <delay>200</delay>
 *  </detector>
 * </root>
 *
 * getChildValue("detector",0,"channel"); will return "C1";
 * and
 * getChildValue("detector",1,"delay"); will return "200";
 *
 * The idea is to then cast the string into the type desired for the child
 */

std::string XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber, std::string _childNode){

      XMLCh* bufferNode = XMLString::transcode(_baseNode.c_str());
      DOMNodeList* list = myDOMParser->getDocument()->getElementsByTagName(bufferNode);
      XMLString::release(&bufferNode);

      if (getBaseNodeCount(_baseNode) <= _baseNumber)
        return "";

      DOMElement* parent = dynamic_cast<DOMElement*>(list->item(_baseNumber));
      DOMElement* child = dynamic_cast<DOMElement*>(parent->getElementsByTagName(XMLString::transcode(_childNode.c_str()))->item(0));
      std::string value;
      if (child) {
        char* bufferChild = XMLString::transcode(child->getTextContent());
        value = bufferChild;
        XMLString::release(&bufferChild);
      } else {
        //Not found? then empty
        value = "";
      }
      return value;
}

/**
 * @brief This methods perform the same as XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber, std::string _childNod) but returns the output into a string passed as argument
 * @param _baseNode the name of the baseNode
 * @param _baseNumber the number of the baseNode to choose
 * @param _childNode the name of the childNode to choose
 * @param _retVal the string containing the result
 */
void  XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber,
    std::string _childNode, std::string &_retVal)
{
  _retVal = getChildValue(_baseNode, _baseNumber, _childNode);
}

/**
 * @brief This methods perform the same as XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber, std::string _childNod) but cast the output into a boolean passed as argument
 * @param _baseNode the name of the baseNode
 * @param _baseNumber the number of the baseNode to choose
 * @param _childNode the name of the childNode to choose
 * @param _retVal the boolean containing the result
 */
void  XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber,
    std::string _childNode, bool& _retVal)
{
    _retVal = (getChildValue(_baseNode, _baseNumber, _childNode)!= "false");
}

/**
 * @brief This methods perform the same as XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber, std::string _childNod) but cast the output into a double passed as argument
 * @param _baseNode the name of the baseNode
 * @param _baseNumber the number of the baseNode to choose
 * @param _childNode the name of the childNode to choose
 * @param _retVal the double containing the result
 */
void  XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber,
    std::string _childNode, double& _retVal)
{
    std::stringstream buffer;
    buffer << getChildValue(_baseNode, _baseNumber, _childNode);
    buffer >> _retVal;
}

/**
 * @brief This methods perform the same as XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber, std::string _childNod) but cast the output into a integer passed as argument
 * @param _baseNode the name of the baseNode
 * @param _baseNumber the number of the baseNode to choose
 * @param _childNode the name of the childNode to choose
 * @param _retVal the double containing the result
 */
void  XMLSettingsReader::getChildValue(std::string _baseNode, int _baseNumber,
    std::string _childNode, int& _retVal)
{
    std::stringstream buffer;
    buffer << getChildValue(_baseNode, _baseNumber, _childNode);
    buffer >> _retVal;
}




