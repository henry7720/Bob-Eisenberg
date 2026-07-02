import re
from docx import Document
from docx.oxml.ns import qn

def update_cv_links(docx_path, output_path):
    doc = Document(docx_path)
    
    # FIX: Changed (.*) to ([^\s)\]]+)
    # This forces the capture group to match the filename but STOP if it hits:
    # a space (\s), a closing parenthesis (\)), or a closing bracket (\])
    reprint_pattern = re.compile(
        r'((?:https?:\/\/|ftp:\/\/)?(?:ftp\.rush\.edu)?\/users\/molebio\/Bob_Eisenberg\/Reprints\/)([^\s)\]]+)', 
        re.IGNORECASE
    )
    
    new_base_url = "https://bob-eisenberg.com/Reprints/"
    hyperlink_updates = 0
    text_updates = 0

    # --- Step 1: Update Underlying Hyperlink Destinations (.rels) ---
    rels = doc.part.rels
    for rel_id, rel in rels.items():
        if "hyperlink" in rel.reltype:
            match = reprint_pattern.search(rel.target_ref)
            if match:
                # group(2) is now our safe filename fragment
                file_fragment = match.group(2)
                rel._target = f"{new_base_url}{file_fragment}"
                hyperlink_updates += 1

    # --- Step 2: Update Raw Text and Visible Link Labels ---
    for paragraph in doc.paragraphs:
        for child in paragraph._element.iter():
            if child.tag.endswith('t') and child.text:
                # We use re.sub here to cleanly replace ONLY the matched URL string
                # leaving everything else in the text run (like ' [PDF]') completely untouched
                new_text, count = reprint_pattern.subn(rf"{new_base_url}\2", child.text)
                if count > 0:
                    child.text = new_text
                    text_updates += count

    # --- Step 3: Check inside Tables ---
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    for child in paragraph._element.iter():
                        if child.tag.endswith('t') and child.text:
                            new_text, count = reprint_pattern.subn(rf"{new_base_url}\2", child.text)
                            if count > 0:
                                child.text = new_text
                                text_updates += count

    doc.save(output_path)
    print("--- Link Refactoring Metrics ---")
    print(f"✔ Modified Hyperlinks (Targets): {hyperlink_updates}")
    print(f"✔ Modified Plain Text References: {text_updates}")
    print(f"✔ Successfully compiled document: '{output_path}'")

update_cv_links("../newest-CV/Bob_Eisenberg_CV_2026-02-19-2.docx", "../newest-CV/Bob_Eisenberg_CV_2026-02-19.docx")
