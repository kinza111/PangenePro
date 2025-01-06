import os
import xml.etree.ElementTree as ET

# Directory containing all XML output files
xml_parent_directory = '/PangenePro/IPS_out'

# Output directory for domain analysis results
output_directory = 'DomAnals'
os.makedirs(output_directory, exist_ok=True)

# Open the output file for writing
for subdir in os.listdir(xml_parent_directory):
    if subdir.startswith('out'):
        xml_directory = os.path.join(xml_parent_directory, subdir)
        output_file_path = os.path.join(output_directory, f'protein_domains_all_{subdir[3:]}.txt')
        
        with open(output_file_path, 'w') as output_file:
            # Write header
            output_file.write('Protein ID\tInterPro Domains\tDescription\n')

            # Iterate over all files in the XML directory
            for filename in os.listdir(xml_directory):
                # Check if file is XML
                if filename.endswith('.xml'):
                    xml_file = os.path.join(xml_directory, filename)

                    # Parse the XML file
                    tree = ET.parse(xml_file)
                    root = tree.getroot()

                    # Iterate over protein elements
                    for protein in root.findall('.//{https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas}protein'):
                        # Get protein ID
                        protein_id = protein.find('.//{https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas}xref').attrib['id']

                        # Initialize lists to store domain information
                        domain_accessions = []
                        domain_descriptions = []

                        # Get InterPro domains
                        for entry in protein.findall('.//{https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas}entry'):
                            if entry.attrib['type'] == 'DOMAIN':
                                domain_accessions.append(entry.attrib['ac'])
                                domain_descriptions.append(entry.attrib['desc'])

                        # Write protein ID, InterPro domains, and descriptions to the output file
                        for accession, description in zip(domain_accessions, domain_descriptions):
                            output_file.write(f"{protein_id}\t{accession}\t{description}\n")

print("Extraction complete. Please check the 'DomAnals' directory for the results.")
