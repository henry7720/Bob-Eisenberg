import re
from docx import Document
from docx.oxml.ns import qn

def update_cv_links(docx_path, output_path):
    # Load the Word document
    doc = Document(docx_path)
    
    # Updated regex to explicitly catch the '/users/' path segment and match both http/https/ftp protocols
    reprint_pattern = re.compile(
        r'(?:https?:\/\/|ftp:\/\/)?(?:ftp\.rush\.edu)?\/users\/molebio\/Bob_Eisenberg\/Reprints\/(.*)', 
        re.IGNORECASE
    )
    
    new_base_url = "https://bob-eisenberg.com/Reprints/"
    hyperlink_updates = 0
    text_updates = 0

# --- Step 1: Update Underlying Hyperlink Destinations (.rels) ---
    rels = doc.part.rels
    for rel_id, rel in rels.items():
        if "hyperlink" in rel.reltype:
            # Read the target link using the public attribute
            match = reprint_pattern.search(rel.target_ref)
            if match:
                file_fragment = match.group(1)
                # Write to the protected internal target attribute to bypass the read-only restriction
                rel._target = f"{new_base_url}{file_fragment}"
                hyperlink_updates += 1

    # --- Step 2: Update Raw Text and Visible Link Labels ---
    for paragraph in doc.paragraphs:
        for child in paragraph._element.iter():
            if child.tag.endswith('t') and child.text:
                match = reprint_pattern.search(child.text)
                if match:
                    file_fragment = match.group(1)
                    child.text = f"{new_base_url}{file_fragment}"
                    text_updates += 1

    # --- Step 3: Check inside Tables ---
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    for child in paragraph._element.iter():
                        if child.tag.endswith('t') and child.text:
                            match = reprint_pattern.search(child.text)
                            if match:
                                file_fragment = match.group(1)
                                child.text = f"{new_base_url}{file_fragment}"
                                text_updates += 1

    # Save changes to the newly targeted document path
    doc.save(output_path)
    print("--- Link Refactoring Metrics ---")
    print(f"✔ Modified Hyperlinks (Targets): {hyperlink_updates}")
    print(f"✔ Modified Plain Text References: {text_updates}")
    print(f"✔ Successfully compiled document: '{output_path}'")

# Execute conversion using your exact filename parameters
update_cv_links("../newest-CV/__CV_January_13_2025_(JT).docx", "../newest-CV/Updated___CV_January_13_2025_(JT).docx")
