"""Catalog warnings must not silently erase coverage or corrupt period values."""
import importlib.util
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location("prior_review", Path(__file__).parents[1] / "scripts" / "review_qlp_catalogs.py")
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


def votable(fields, row="", info=""):
    return (f'<?xml version="1.0"?><VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">'
            f'<RESOURCE>{info}<TABLE name="catalog">{fields}'
            f'<DATA><TABLEDATA>{row}</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>').encode()


def test_reviewed_generated_xml_ids_preserve_column_names_and_values():
    body = votable('<FIELD name="u_M-m/D" datatype="double" unit="%"/>'
                   '<FIELD name="E(B-V)" datatype="double"/>'
                   '<FIELD name="2M" datatype="int"/>',
                   '<TR><TD>12.5</TD><TD>0.02</TD><TD>7</TD></TR>')
    parsed, recorded = review.parse_vizier(body)
    table = next(parsed.iter_tables())
    assert [field.name for field in table.fields] == ["u_M-m/D", "E(B-V)", "2M"]
    assert tuple(table.array[0]) == (12.5, 0.02, 7)
    assert [item["code"] for item in recorded].count("W03") == 3
    assert [item["code"] for item in recorded].count("W50") == 1


def test_unreviewed_unit_warning_is_a_failure():
    with pytest.raises(ValueError, match="Unreviewed VOTable verification warning"):
        review.parse_vizier(votable('<FIELD name="period" datatype="double" unit="bad!unit"/>'))


@pytest.mark.parametrize("value", ["ERROR", "OVERFLOW"])
def test_service_errors_and_truncation_fail(value):
    with pytest.raises(ValueError, match="VizieR errors/overflow"):
        review.parse_vizier(votable('<FIELD name="period" datatype="double" unit="d"/>',
                                   info=f'<INFO name="QUERY_STATUS" value="{value}"/>'))
