#!/usr/bin/env python3
"""Validate an XML file against a schema.

Written by William E Fondrie, 2021-05-19
"""
import sys
import logging
from argparse import ArgumentParser

try:
    from lxml import etree
except ImportError as err:
    raise ImportError(
        "Please install the lxml package with 'pip install lxml' or 'conda "
        "install lxml'."
    )

LOGGER = logging.getLogger(__name__)

# The URLs for the PepXML and mzIdentML schemas. Update these when new versions
# are released.
PEPXML_URL = "http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v122.xsd"
MZID_URL = "https://www.psidev.info/sites/default/files/2017-06/mzIdentML1.1.0.xsd"


def parse_args():
    """Parse command line arguments

    Returns
    -------
    argparse.ArgumentParser
        An argument parser with the values of each argument as an attribute.
    """
    parser = ArgumentParser(
        description=(
            "Validate an XML file against an XML schema. All errors that are "
            "found are written to stdout."
        )
    )
    parser.add_argument("xml_path", help="The XML file to validate.")
    parser.add_argument(
        "-s",
        "--xsd_path",
        default="pepxml",
        help=(
            "The XML schema to validate against. This can be a file, URL, or "
            "one of {pepxml, mzid}. The latter will automatically grab the "
            "PepXML and mzIdentML format schemas from their respective URLs."
        ),
    )
    return parser.parse_args()


def find_xsd(xsd_path):
    """Find the PepXML or mzIdentML file schemas if requested.

    Parameters
    ----------
    xsd_path : str or Path
        "pxpxml", "mzid", or the XML schema to validate against.

    Returns
    -------
    str
        The path to the XML schema, local or otherwise.
    """
    if xsd_path.lower() == "pepxml":
        return PEPXML_URL
    if xsd_path.lower() == "mzid" or xsd_path.lower() == "mzidentml":
        return MZID_URL

    return xsd_path


def validate(xml_path, xsd_path):
    """Validate an XML file using the XSD schema.

    Parameters
    ----------
    xml_path : str or Path
        The XML file to validate.
    xsd_path : str of Path
        The XML schema to validate against.
    """
    LOGGER.info("Parsing the XSD file: %s", xsd_path)
    xsd = etree.XMLSchema(etree.parse(xsd_path))

    LOGGER.info("Parsing the XML file: %s", xml_path)
    xml = etree.parse(xml_path)

    LOGGER.info("Validating...")
    try:
        xsd.assertValid(xml)
    except etree.DocumentInvalid as err:
        sys.stdout.write(str(err.error_log))
        raise etree.DocumentInvalid("The XML file failed to validate.")


def main():
    """The main function"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    args = parse_args()
    xsd_path = find_xsd(args.xsd_path)
    validate(args.xml_path, xsd_path)
    LOGGER.info("Successfully validated the XML file.")


if __name__ == "__main__":
    main()
