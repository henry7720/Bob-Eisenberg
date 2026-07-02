import re
from docx import Document

NEW_BASE_URL = "https://bob-eisenberg.com/Reprints/"

REPRINT_PATTERN = re.compile(
    r'(?:https?://|ftp://)?'
    r'(?:ftp\.rush\.edu)?'
    r'/?users/molebio/Bob_Eisenberg/Reprints/'
    r'([^\s)\]]+)',
    re.IGNORECASE,
)


def replacement(match):
    """Return the new URL preserving only the filename/path fragment."""
    return NEW_BASE_URL + match.group(1)


def update_text_in_paragraphs(paragraphs):
    """Update visible URLs in a collection of paragraphs."""
    updates = 0

    for paragraph in paragraphs:
        for node in paragraph._element.iter():
            if node.tag.endswith("t") and node.text:
                new_text, count = REPRINT_PATTERN.subn(replacement, node.text)
                if count:
                    node.text = new_text
                    updates += count

    return updates


def update_tables(tables):
    """Recursively update all tables."""
    updates = 0

    for table in tables:
        for row in table.rows:
            for cell in row.cells:
                updates += update_text_in_paragraphs(cell.paragraphs)

                # Handle nested tables
                updates += update_tables(cell.tables)

    return updates


def update_cv_links(docx_path, output_path):
    doc = Document(docx_path)

    hyperlink_updates = 0
    text_updates = 0

    # ------------------------------------------------------------------
    # Update hyperlink destinations (.rels)
    # ------------------------------------------------------------------
    for rel in doc.part.rels.values():
        if rel.reltype.endswith("/hyperlink"):
            match = REPRINT_PATTERN.search(rel.target_ref)
            if match:
                rel._target = replacement(match)
                hyperlink_updates += 1

    # ------------------------------------------------------------------
    # Main document
    # ------------------------------------------------------------------
    text_updates += update_text_in_paragraphs(doc.paragraphs)
    text_updates += update_tables(doc.tables)

    # ------------------------------------------------------------------
    # Headers and footers
    # ------------------------------------------------------------------
    for section in doc.sections:
        header = section.header
        footer = section.footer

        text_updates += update_text_in_paragraphs(header.paragraphs)
        text_updates += update_tables(header.tables)

        text_updates += update_text_in_paragraphs(footer.paragraphs)
        text_updates += update_tables(footer.tables)

    doc.save(output_path)

    print("\n--- Link Refactoring Metrics ---")
    print(f"✔ Modified hyperlink targets : {hyperlink_updates}")
    print(f"✔ Modified visible URLs      : {text_updates}")
    print(f"✔ Saved updated document to: {output_path}")


if __name__ == "__main__":
    update_cv_links(
        "../newest-CV/Bob_Eisenberg_CV_2026-02-19-2.docx",
        "../newest-CV/Bob_Eisenberg_CV_2026-02-19-3.docx",
    )
