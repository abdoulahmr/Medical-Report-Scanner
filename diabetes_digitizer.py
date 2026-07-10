#!/usr/bin/env python3

import os
import re
import csv
import json
import glob
import queue
import string
import hashlib
import difflib
import threading
import unicodedata
from datetime import datetime

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import cv2
from PIL import Image, ImageTk

try:
    from paddleocr import PaddleOCR
except ImportError:
    PaddleOCR = None


IMAGE_EXTENSIONS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff", ".webp"}

KNOWN_FIELDS = {
    "glycemie a jeun": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "glycemie": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "glycemie post prandiale": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "glycemie postprandiale": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "hba1c": {"unit": "%", "min": 3.0, "max": 18.0},
    "hemoglobine glyquee": {"unit": "%", "min": 3.0, "max": 18.0},
    "poids": {"unit": "kg", "min": 25.0, "max": 220.0},
    "taille": {"unit": "cm", "min": 100.0, "max": 210.0},
    "imc": {"unit": "kg/m2", "min": 12.0, "max": 60.0},
    "tension arterielle systolique": {"unit": "mmHg", "min": 70.0, "max": 260.0},
    "tension arterielle diastolique": {"unit": "mmHg", "min": 40.0, "max": 150.0},
    "creatinine": {"unit": "mg/L", "min": 2.0, "max": 100.0},
    "cholesterol total": {"unit": "g/L", "min": 0.5, "max": 5.0},
    "triglycerides": {"unit": "g/L", "min": 0.2, "max": 10.0},
}

LABEL_SUGGESTIONS = [
    "HEMATIES", "Hémoglobine", "Hematocrite","V.G.M", "T.C.M.H", "C.C.M.H.",
    "I.D.R(R.D.W)", "LEUCOCYTES","Polynucléaires Neutrophiles %",
    "Polynucléaires Neutrophiles mm3", "Lymphocytes %", "Lymphocytes mm3", 
    "Monocytes %", "Monocytes mm3", "Eosinophiles %", "Eosinophiles mm3",
    "Polynucléaires Basophiles %", "Polynucléaires Basophiles mm3",
    "Glycime a jeun g/l", "Glycime a jeun mmol/l", "Uree sanguine", 
    "Creatinine sanguine", "Transaminases SGOT", "Transaminases SGPT", 
    "Hemoglobine glyquee HbA1c","PLAQUETTES", "MPV", "PDW", "P-LCR", 
    "Volume Urinaire des 24 heures", "CHOLESTEROLTOTAL g/l", "CHOLESTEROLTOTAL mmol/l", 
    "TRIGLYCERIDES g/l", "TRIGLYCERIDES mmol/l","CHOLESTEROL HDL", "CHOLESTEROL LDL", "RAPPORT Cholesterol T/HDL Chol",
    "MICROALBUMINURIE mg/l", "MICROALBUMINURIE mg/24h", "T.S.H. ultra-sensible"
]

# Different labs use different acronyms/abbreviations for the same test
# (e.g. "Gly" or "GAJ" for "Glycime a jeun g/l"). This is a STARTER set of
# common French lab abbreviations - it's deliberately not exhaustive.
# Anything not covered here gets learned automatically the first time you
# correct a label (see DataStore.learn_alias), or can be added by hand via
# the "Manage Synonyms..." button.
BUILTIN_ALIASES = {
    "HEMATIES": ["GR", "Hematies", "RBC"],
    "Hémoglobine": ["Hb", "HGB", "Hb g/dl"],
    "Hematocrite": ["Hte", "HCT", "Ht"],
    "V.G.M": ["VGM", "MCV"],
    "T.C.M.H": ["TCMH", "MCH"],
    "C.C.M.H.": ["CCMH", "MCHC"],
    "I.D.R(R.D.W)": ["IDR", "RDW"],
    "LEUCOCYTES": ["GB", "WBC", "Leuco"],
    "Polynucléaires Neutrophiles %": ["PNN %", "Neutro %", "Neut %"],
    "Polynucléaires Neutrophiles mm3": ["PNN mm3", "Neutro mm3", "Neut mm3"],
    "Lymphocytes %": ["Lympho %", "Lym %"],
    "Lymphocytes mm3": ["Lympho mm3", "Lym mm3"],
    "Monocytes %": ["Mono %"],
    "Monocytes mm3": ["Mono mm3"],
    "Eosinophiles %": ["Eosino %", "Eo %"],
    "Eosinophiles mm3": ["Eosino mm3", "Eo mm3"],
    "Polynucléaires Basophiles %": ["PNB %", "Baso %"],
    "Polynucléaires Basophiles mm3": ["PNB mm3", "Baso mm3"],
    "Glycime a jeun g/l": ["Gly", "Glyc", "GAJ", "Glycemie"],
    "Glycime a jeun mmol/l": ["Gly mmol", "GAJ mmol"],
    "Uree sanguine": ["Uree", "BUN"],
    "Creatinine sanguine": ["Creat", "Crea", "Creatinine"],
    "Transaminases SGOT": ["SGOT", "ASAT", "TGO"],
    "Transaminases SGPT": ["SGPT", "ALAT", "TGP"],
    "Hemoglobine glyquee HbA1c": ["HbA1c", "Hb glyquee", "Hb A1c"],
    "PLAQUETTES": ["Plaq", "PLT"],
    "CHOLESTEROLTOTAL g/l": ["Chol T", "CT", "Cholesterol total"],
    "TRIGLYCERIDES g/l": ["TG", "Triglycerides"],
    "CHOLESTEROL HDL": ["HDL", "HDL-C", "HDLc"],
    "CHOLESTEROL LDL": ["LDL", "LDL-C", "LDLc"],
    "MICROALBUMINURIE mg/l": ["MAU", "Microalb"],
    "T.S.H. ultra-sensible": ["TSH", "TSHus"],
}

DATE_PATTERNS = [
    re.compile(r"\b(\d{1,2})[/\-.](\d{1,2})[/\-.](\d{2,4})\b"),
]


def normalize_label(label):
    label = label.strip().lower()
    label = unicodedata.normalize("NFKD", label)
    label = "".join(c for c in label if not unicodedata.combining(c))
    label = re.sub(r"\s+", " ", label)
    return label


def try_parse_date(text):
    for pat in DATE_PATTERNS:
        m = pat.search(text)
        if m:
            d, mth, y = m.groups()
            try:
                d, mth, y = int(d), int(mth), int(y)
                if y < 100:
                    y += 2000 if y < 50 else 1900
                dt = datetime(y, mth, d)
                return dt.strftime("%Y-%m-%d")
            except ValueError:
                continue
    return None


# --------------------------------------------------------------------------
# Row reconstruction: PaddleOCR often returns one box per "word cluster",
# so a printed line like
#     HEMATIES ..................... 4,90 10^6/mm3
# comes back as THREE separate boxes (label / value / unit) because the
# dot leaders create a gap the detector doesn't bridge. We reconstruct the
# original line by grouping boxes with similar Y position into rows, then
# sorting each row left-to-right and concatenating the text.
# --------------------------------------------------------------------------

DOT_LEADER_RE = re.compile(r"(?:\.\s*){2,}")
NUMBER_RE = re.compile(r"-?\d+(?:[.,]\d+)?")

_LABEL_LOOKUP_NORM = {}
for _suggestion in LABEL_SUGGESTIONS:
    _LABEL_LOOKUP_NORM[normalize_label(_suggestion)] = _suggestion

_BUILTIN_ALIAS_LOOKUP_NORM = {}
for _canonical, _synonyms in BUILTIN_ALIASES.items():
    for _syn in _synonyms:
        _BUILTIN_ALIAS_LOOKUP_NORM[normalize_label(_syn)] = _canonical


def group_ocr_lines(raw_boxes):
    """Group raw OCR boxes into reconstructed text rows.

    Steps: sort by Y center, cluster boxes whose Y center is close together
    into the same row, sort each row's boxes by X, then concatenate text.
    Returns a list of {"text", "confidence", "bbox"} dicts, one per row,
    sorted top-to-bottom.
    """
    if not raw_boxes:
        return []

    items = []
    for box in raw_boxes:
        xs = [p[0] for p in box["bbox"]]
        ys = [p[1] for p in box["bbox"]]
        items.append({
            "text": box["text"],
            "confidence": box["confidence"],
            "x0": min(xs), "x1": max(xs),
            "y0": min(ys), "y1": max(ys),
            "cy": (min(ys) + max(ys)) / 2.0,
            "h": max(1.0, max(ys) - min(ys)),
        })
    items.sort(key=lambda it: it["cy"])

    rows = []  # each: {"items": [...], "cy_sum": float, "h_sum": float, "count": int}
    for it in items:
        placed = False
        for row in rows:
            row_cy = row["cy_sum"] / row["count"]
            row_h = row["h_sum"] / row["count"]
            tolerance = max(row_h, it["h"]) * 0.7
            if abs(it["cy"] - row_cy) <= tolerance:
                row["items"].append(it)
                row["cy_sum"] += it["cy"]
                row["h_sum"] += it["h"]
                row["count"] += 1
                placed = True
                break
        if not placed:
            rows.append({"items": [it], "cy_sum": it["cy"], "h_sum": it["h"], "count": 1})

    grouped = []
    for row in rows:
        row_items = sorted(row["items"], key=lambda it: it["x0"])
        text = " ".join(it["text"] for it in row_items)
        text = DOT_LEADER_RE.sub(" ", text)
        text = re.sub(r"\s+", " ", text).strip()
        if not text:
            continue
        x0 = min(it["x0"] for it in row_items)
        x1 = max(it["x1"] for it in row_items)
        y0 = min(it["y0"] for it in row_items)
        y1 = max(it["y1"] for it in row_items)
        conf = sum(it["confidence"] for it in row_items) / len(row_items)
        grouped.append({
            "text": text,
            "confidence": conf,
            "bbox": [[x0, y0], [x1, y0], [x1, y1], [x0, y1]],
        })

    grouped.sort(key=lambda g: g["bbox"][0][1])
    return grouped


def parse_label_value_unit(text):
    """Split a reconstructed line like 'HEMATIES 4,90 10^6/mm3' into
    (label, value, unit) using the first standalone number as the split
    point. Returns ("", "", "") pieces where nothing plausible is found."""
    cleaned = DOT_LEADER_RE.sub(" ", text)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    m = NUMBER_RE.search(cleaned)
    if not m:
        return cleaned, "", ""
    label = cleaned[:m.start()].strip(" .:-\u2013\u2014")
    value = m.group(0)
    unit = cleaned[m.end():].strip(" .:-\u2013\u2014")
    return label, value, unit


def match_known_label(raw_label, extra_aliases=None):
    """Resolve an OCR'd label to a canonical LABEL_SUGGESTIONS entry.

    Different labs abbreviate the same test differently (e.g. "Gly" for
    "Glycime a jeun g/l"), and an abbreviation rarely looks similar to the
    full name character-by-character, so fuzzy string matching alone won't
    catch it. Resolution order:
      1. Exact match against a canonical label.
      2. Exact match against a known alias (built-in list, or ones learned
         from your past corrections in this folder via extra_aliases).
      3. Fuzzy match against canonical labels only (handles OCR typos on
         the full name; NOT used for aliases, since short abbreviations
         are too easy to fuzzy-match to the wrong thing).
      4. Otherwise, return the raw OCR text unchanged so it's still
         editable in the review list.
    """
    if not raw_label:
        return raw_label
    norm = normalize_label(raw_label)
    if not norm:
        return raw_label

    if norm in _LABEL_LOOKUP_NORM:
        return _LABEL_LOOKUP_NORM[norm]

    if norm in _BUILTIN_ALIAS_LOOKUP_NORM:
        return _BUILTIN_ALIAS_LOOKUP_NORM[norm]

    if extra_aliases and norm in extra_aliases:
        return extra_aliases[norm]["canonical"]

    best = difflib.get_close_matches(norm, _LABEL_LOOKUP_NORM.keys(), n=1, cutoff=0.7)
    if best:
        return _LABEL_LOOKUP_NORM[best[0]]

    return raw_label


def validate_value(label, value_text):
    key = normalize_label(label)
    spec = KNOWN_FIELDS.get(key)
    if not spec:
        return None
    num_match = re.search(r"-?\d+[.,]?\d*", value_text)
    if not num_match:
        return "expected a numeric value"
    try:
        val = float(num_match.group(0).replace(",", "."))
    except ValueError:
        return None
    if val < spec["min"] or val > spec["max"]:
        return f"unusual value for {label} (typical range {spec['min']}-{spec['max']} {spec['unit']})"
    return None


def file_hash(path):
    st = os.stat(path)
    key = f"{os.path.abspath(path)}|{st.st_size}|{int(st.st_mtime)}"
    return hashlib.sha1(key.encode("utf-8")).hexdigest()


def order_points(pts):
    """Order 4 points as top-left, top-right, bottom-right, bottom-left."""
    rect = np.zeros((4, 2), dtype="float32")
    s = pts.sum(axis=1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]
    diff = np.diff(pts, axis=1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]
    return rect


def find_document_contour(img):
    """Find the largest 4-point contour that plausibly outlines a
    document/page in the image. Returns a (4, 2) float32 array of corner
    points, or None if nothing clearly document-shaped was found."""
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    edged = cv2.Canny(blurred, 50, 150)
    edged = cv2.dilate(edged, np.ones((5, 5), np.uint8), iterations=2)
    edged = cv2.erode(edged, np.ones((5, 5), np.uint8), iterations=1)

    contours, _ = cv2.findContours(edged, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None
    contours = sorted(contours, key=cv2.contourArea, reverse=True)[:10]

    img_area = img.shape[0] * img.shape[1]
    for c in contours:
        peri = cv2.arcLength(c, True)
        approx = cv2.approxPolyDP(c, 0.02 * peri, True)
        if len(approx) == 4 and cv2.contourArea(approx) > 0.2 * img_area:
            return approx.reshape(4, 2).astype("float32")
    return None


def four_point_transform(img, pts):
    """Warp the quadrilateral defined by `pts` into a flat, top-down
    rectangle (standard perspective-correction / "scan" transform)."""
    rect = order_points(pts)
    (tl, tr, br, bl) = rect

    width_a = np.linalg.norm(br - bl)
    width_b = np.linalg.norm(tr - tl)
    max_width = max(int(width_a), int(width_b), 1)

    height_a = np.linalg.norm(tr - br)
    height_b = np.linalg.norm(tl - bl)
    max_height = max(int(height_a), int(height_b), 1)

    dst = np.array([
        [0, 0],
        [max_width - 1, 0],
        [max_width - 1, max_height - 1],
        [0, max_height - 1],
    ], dtype="float32")

    m = cv2.getPerspectiveTransform(rect, dst)
    return cv2.warpPerspective(img, m, (max_width, max_height))


def compute_skew_angle(img, max_correction=15.0):
    """Estimate residual (fine) skew angle of a mostly-rectangular document
    from the spread of its dark (text/ink) pixels. Returns degrees; angles
    larger than `max_correction` are treated as a detection failure (e.g.
    a mostly-blank page) and ignored, returning 0.0."""
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    gray = cv2.bitwise_not(gray)
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]

    coords = np.column_stack(np.where(thresh > 0))
    if coords.shape[0] < 20:
        return 0.0

    angle = cv2.minAreaRect(coords)[-1]
    if angle < -45:
        angle = -(90 + angle)
    else:
        angle = -angle

    if abs(angle) > max_correction:
        return 0.0
    return angle


def remove_shadows(img):
    """Flatten uneven lighting/shadows on a scanned document.

    For each color channel, estimates the background illumination via
    dilation + median blur, then normalizes the image against that
    background so shadows and gradients are removed while text stays dark.
    """
    planes = cv2.split(img)
    result_planes = []
    for plane in planes:
        dilated = cv2.dilate(plane, np.ones((7, 7), np.uint8))
        bg = cv2.medianBlur(dilated, 21)
        diff = 255 - cv2.absdiff(plane, bg)
        norm = cv2.normalize(diff, None, alpha=0, beta=255,
                              norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_8UC1)
        result_planes.append(norm)
    return cv2.merge(result_planes)


class ImagePipeline:

    def __init__(self, path):
        self.path = path
        self.original = cv2.imread(path, cv2.IMREAD_COLOR)
        if self.original is None:
            raise ValueError(f"Could not read image: {path}")
        self.rotation_90 = 0
        self.fine_angle = 0.0
        self.crop_rect = None
        self.denoise = False
        self.contrast = False
        # Set by auto_scan(): a perspective-corrected, deskewed version of
        # `original`. When present, render() builds on top of this instead
        # of the raw original. None means "no scan run yet".
        self.scanned_base = None
        self.shadow_removal = False

    def reset(self):
        self.rotation_90 = 0
        self.fine_angle = 0.0
        self.crop_rect = None
        self.denoise = False
        self.contrast = False
        self.scanned_base = None
        self.shadow_removal = False

    def auto_scan(self):
        """One-click document scan: detect the page edges, perspective-
        correct it flat, auto-rotate to portrait, deskew any residual tilt,
        remove shadows, and enable contrast enhancement.

        Replaces `scanned_base` (manual crop/rotate/angle adjustments are
        cleared since the scan already reframes the page); denoise and the
        shadow-removal/contrast toggles remain available afterwards.

        Returns True if a document boundary was detected, False if the
        scan fell back to using the full original frame.
        """
        self.crop_rect = None
        self.rotation_90 = 0
        self.fine_angle = 0.0

        working = self.original.copy()

        quad = find_document_contour(working)
        found = quad is not None
        scanned = four_point_transform(working, quad) if found else working

        # Auto-rotate: FNS/lab forms are portrait; a landscape result after
        # perspective correction almost always means the page was
        # photographed sideways.
        h, w = scanned.shape[:2]
        if w > h * 1.15:
            scanned = cv2.rotate(scanned, cv2.ROTATE_90_CLOCKWISE)

        skew = compute_skew_angle(scanned)
        if abs(skew) > 0.3:
            hh, ww = scanned.shape[:2]
            center = (ww / 2, hh / 2)
            m = cv2.getRotationMatrix2D(center, skew, 1.0)
            scanned = cv2.warpAffine(scanned, m, (ww, hh),
                                      flags=cv2.INTER_LINEAR,
                                      borderMode=cv2.BORDER_REPLICATE)

        self.scanned_base = scanned
        self.shadow_removal = True
        self.contrast = True
        return found

    def render(self):
        img = (self.scanned_base if self.scanned_base is not None else self.original).copy()

        if self.crop_rect:
            x1, y1, x2, y2 = self.crop_rect
            x1, x2 = sorted((max(0, x1), max(0, x2)))
            y1, y2 = sorted((max(0, y1), max(0, y2)))
            h, w = img.shape[:2]
            x2, y2 = min(w, x2), min(h, y2)
            if x2 - x1 > 5 and y2 - y1 > 5:
                img = img[y1:y2, x1:x2]

        if self.rotation_90 == 90:
            img = cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE)
        elif self.rotation_90 == 180:
            img = cv2.rotate(img, cv2.ROTATE_180)
        elif self.rotation_90 == 270:
            img = cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)

        if abs(self.fine_angle) > 0.01:
            h, w = img.shape[:2]
            center = (w / 2, h / 2)
            m = cv2.getRotationMatrix2D(center, self.fine_angle, 1.0)
            img = cv2.warpAffine(img, m, (w, h),
                                  flags=cv2.INTER_LINEAR,
                                  borderMode=cv2.BORDER_REPLICATE)

        if self.shadow_removal:
            if len(img.shape) == 2:
                img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            elif img.shape[2] == 4:
                img = cv2.cvtColor(img, cv2.COLOR_BGRA2BGR)
            img = remove_shadows(img)

        if self.denoise:
            if len(img.shape) == 2:
                img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            elif img.shape[2] == 4:
                img = cv2.cvtColor(img, cv2.COLOR_BGRA2BGR)
            img = cv2.fastNlMeansDenoisingColored(
                src=img,
                dst=None,
                h=7,
                hColor=7,
                templateWindowSize=7,
                searchWindowSize=21
            )

        if self.contrast:
            lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
            l, a, b = cv2.split(lab)
            clahe = cv2.createCLAHE(clipLimit=2.5, tileGridSize=(8, 8))
            l = clahe.apply(l)
            lab = cv2.merge((l, a, b))
            img = cv2.cvtColor(lab, cv2.COLOR_LAB2BGR)

        return img


class OCREngine:

    def __init__(self):
        self._ocr = None
        self._lock = threading.Lock()

    @staticmethod
    def _paddleocr_cache_dir():
        return os.path.join(os.path.expanduser("~"), ".paddleocr")

    def _init_paddleocr(self):
        return PaddleOCR(use_angle_cls=True, lang="fr", show_log=False)

    def ensure_loaded(self):
        with self._lock:
            if self._ocr is None:
                if PaddleOCR is None:
                    raise RuntimeError(
                        "PaddleOCR is not installed. Run:\n"
                        "  pip install paddleocr paddlepaddle"
                    )
                try:
                    self._ocr = self._init_paddleocr()
                except Exception as e:
                    msg = str(e).lower()
                    corruption_hints = ("unexpected end of data", "tarfile",
                                         "eof", "not a gzip file", "truncat")
                    cache_dir = self._paddleocr_cache_dir()
                    if any(h in msg for h in corruption_hints) and os.path.isdir(cache_dir):
                        import shutil
                        shutil.rmtree(cache_dir, ignore_errors=True)
                        try:
                            self._ocr = self._init_paddleocr()
                        except Exception as e2:
                            raise RuntimeError(
                                "PaddleOCR model files looked corrupted, so the "
                                f"cache at {cache_dir} was cleared and a fresh "
                                "download was attempted, but it still failed:\n"
                                f"{e2}\n\n"
                                "Check your internet connection and try Run OCR again."
                            ) from e2
                    else:
                        raise
        return self._ocr

    def run(self, bgr_image):
        ocr = self.ensure_loaded()
        result = ocr.ocr(bgr_image, cls=True)
        lines = []
        if not result:
            return lines
        page = result[0] if isinstance(result[0], list) or result[0] is None else result
        page = page or []
        for entry in page:
            try:
                bbox, (text, conf) = entry
            except (ValueError, TypeError):
                continue
            lines.append({
                "bbox": [[float(x), float(y)] for x, y in bbox],
                "text": text,
                "confidence": float(conf),
            })
        return lines


class DataStore:
    def __init__(self, folder):
        self.folder = folder
        self.records_path = os.path.join(folder, "records.json")
        self.manifest_path = os.path.join(folder, "manifest.json")
        self.csv_path = os.path.join(folder, "records.csv")
        self.patients_dir = os.path.join(folder, "patients")
        self.aliases_path = os.path.join(folder, "label_aliases.json")
        self.data = {"patients": {}}
        self.manifest = {"processed": {}}
        self.aliases = {}   # normalized_alias -> {"canonical": str, "example": str}
        self._load()

    def _load(self):
        if os.path.exists(self.records_path):
            with open(self.records_path, "r", encoding="utf-8") as f:
                self.data = json.load(f)
        if os.path.exists(self.manifest_path):
            with open(self.manifest_path, "r", encoding="utf-8") as f:
                self.manifest = json.load(f)
        if os.path.exists(self.aliases_path):
            with open(self.aliases_path, "r", encoding="utf-8") as f:
                self.aliases = json.load(f)
        self._migrate()

    def _migrate(self):
        for p in self.data.get("patients", {}).values():
            for visit in p.get("visits", []):
                if "images" not in visit:
                    old = visit.pop("source_image", None)
                    visit["images"] = [old] if old else []

    def learn_alias(self, raw_text, canonical):
        """Remember that `raw_text` (as typed/OCR'd) means `canonical`, so
        future images with the same abbreviation auto-fill correctly.
        Returns True if this added/changed a mapping."""
        norm = normalize_label(raw_text)
        if not norm or normalize_label(canonical) == norm:
            return False
        existing = self.aliases.get(norm)
        if existing and existing.get("canonical") == canonical:
            return False
        self.aliases[norm] = {"canonical": canonical, "example": raw_text}
        return True

    def remove_alias(self, norm_key):
        self.aliases.pop(norm_key, None)

    def save_aliases(self):
        os.makedirs(self.folder, exist_ok=True)
        with open(self.aliases_path, "w", encoding="utf-8") as f:
            json.dump(self.aliases, f, ensure_ascii=False, indent=2)

    def save(self):
        os.makedirs(self.folder, exist_ok=True)
        os.makedirs(self.patients_dir, exist_ok=True)
        with open(self.records_path, "w", encoding="utf-8") as f:
            json.dump(self.data, f, ensure_ascii=False, indent=2)
        with open(self.manifest_path, "w", encoding="utf-8") as f:
            json.dump(self.manifest, f, ensure_ascii=False, indent=2)
        self.save_aliases()
        self._export_csv()
        self._export_patient_json()

    def is_processed(self, path):
        return file_hash(path) in self.manifest["processed"]

    def find_by_image(self, path):
        entry = self.manifest["processed"].get(file_hash(path))
        if not entry:
            return None, None
        return entry["patient_id"], entry["visit_index"]

    def patient_ids(self):
        return sorted(self.data["patients"].keys())

    def patient_name(self, pid):
        p = self.data["patients"].get(pid)
        return p["name"] if p else ""

    def ensure_patient(self, pid, name):
        if pid not in self.data["patients"]:
            self.data["patients"][pid] = {"name": name, "visits": []}
        elif name:
            self.data["patients"][pid]["name"] = name

    def visits(self, pid):
        p = self.data["patients"].get(pid)
        return p["visits"] if p else []

    def record_options(self, pid):
        options = []
        for idx, visit in enumerate(self.visits(pid)):
            date = visit.get("visit_date") or "no date"
            n_images = len(visit.get("images", []))
            n_fields = len(visit.get("fields", []))
            options.append((idx, f"{date}  -  {n_images} image(s), {n_fields} field(s)"))
        options.sort(key=lambda t: self.visits(pid)[t[0]].get("visit_date") or "")
        return options

    def get_visit(self, pid, visit_index):
        visits = self.visits(pid)
        if visit_index is None or not (0 <= visit_index < len(visits)):
            return None
        return visits[visit_index]

    def save_or_update_visit(self, pid, name, visit_index, image_path, visit_date, fields):
        self.ensure_patient(pid, name)
        visits = self.data["patients"][pid]["visits"]
        abs_path = os.path.abspath(image_path)

        if visit_index is None:
            visit = {
                "images": [abs_path],
                "visit_date": visit_date,
                "processed_at": datetime.now().isoformat(timespec="seconds"),
                "fields": fields,
            }
            visits.append(visit)
            visit_index = len(visits) - 1
        else:
            visit = visits[visit_index]
            if abs_path not in visit["images"]:
                visit["images"].append(abs_path)
            if visit_date:
                visit["visit_date"] = visit_date
            visit["fields"] = fields
            visit["processed_at"] = datetime.now().isoformat(timespec="seconds")

        h = file_hash(image_path)
        self.manifest["processed"][h] = {
            "path": abs_path,
            "patient_id": pid,
            "visit_index": visit_index,
        }
        return visit_index

    def _export_csv(self):
        with open(self.csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["patient_id", "patient_name", "visit_date",
                        "source_images", "label", "value", "confidence",
                        "unit", "flag"])
            for pid, p in sorted(self.data["patients"].items()):
                for visit in p["visits"]:
                    images_joined = "; ".join(os.path.basename(i) for i in visit.get("images", []))
                    for field in visit["fields"]:
                        w.writerow([
                            pid, p["name"], visit["visit_date"],
                            images_joined,
                            field.get("label", ""),
                            field.get("value", ""),
                            field.get("confidence", ""),
                            field.get("unit", ""),
                            field.get("flag", ""),
                        ])

    def _export_patient_json(self):
        for pid, p in self.data["patients"].items():
            path = os.path.join(self.patients_dir, f"{pid}.json")
            visits = sorted(p["visits"], key=lambda v: v.get("visit_date") or "")
            payload = {
                "patient_id": pid,
                "name": p.get("name", ""),
                "visits": [
                    {
                        "visit_date": v.get("visit_date"),
                        "images": [os.path.basename(i) for i in v.get("images", [])],
                        "processed_at": v.get("processed_at"),
                        "fields": v.get("fields", []),
                    }
                    for v in visits
                ],
            }
            with open(path, "w", encoding="utf-8") as f:
                json.dump(payload, f, ensure_ascii=False, indent=2)


class App(tk.Tk):
    CONF_GOOD = 0.85
    CONF_OK = 0.60

    def __init__(self):
        super().__init__()
        self.title("Diabetes Records Digitizer")
        self.geometry("1300x850")

        self.ocr_engine = OCREngine()
        self.store = None
        self.image_paths = []
        self.current_index = -1
        self.pipeline = None
        self.display_scale = 1.0
        self.fit_scale = 1.0
        self.zoom_factor = 1.0
        self.disp_w = 1
        self.disp_h = 1
        self.tk_image = None
        self.ocr_lines = []
        self.grouped_lines = []
        self.field_rows = []
        self.crop_start = None
        self.crop_rect_id = None
        self.ocr_queue = queue.Queue()
        self.current_visit_index = None
        self.record_index_map = []
        self.pending_preselect = None

        # Shared by both the menu bar (checkbuttons/sliders reference them)
        # and the controls row, so they're created before either is built.
        self.angle_var = tk.DoubleVar(value=0.0)
        self.denoise_var = tk.BooleanVar(value=False)
        self.contrast_var = tk.BooleanVar(value=False)
        self.shadow_var = tk.BooleanVar(value=False)
        self.grid_var = tk.BooleanVar(value=False)

        self._build_menu_bar()
        self._build_ui()

    def _build_menu_bar(self):
        menubar = tk.Menu(self, tearoff=0)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Select Folder...", command=self.select_folder,
                               accelerator="Ctrl+O")
        file_menu.add_command(label="Manage Synonyms...", command=self.open_alias_manager)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.destroy, accelerator="Ctrl+Q")
        menubar.add_cascade(label="File", menu=file_menu)

        image_menu = tk.Menu(menubar, tearoff=0)
        image_menu.add_command(label="Scan", command=self.auto_scan_current,
                                accelerator="Ctrl+K")
        image_menu.add_separator()
        image_menu.add_command(label="Rotate 90° Clockwise", command=lambda: self.rotate(90))
        image_menu.add_command(label="Rotate 90° Counter-clockwise", command=lambda: self.rotate(-90))
        image_menu.add_separator()
        image_menu.add_checkbutton(label="Denoise", variable=self.denoise_var,
                                    command=self.refresh_display)
        image_menu.add_checkbutton(label="Contrast (CLAHE)", variable=self.contrast_var,
                                    command=self.refresh_display)
        image_menu.add_checkbutton(label="Remove Shadows", variable=self.shadow_var,
                                    command=self.refresh_display)
        image_menu.add_checkbutton(label="Show Grid", variable=self.grid_var,
                                    command=self.refresh_display)
        image_menu.add_separator()
        image_menu.add_command(label="Clear Crop", command=self.clear_crop)
        image_menu.add_command(label="Reset All Adjustments", command=self.reset_pipeline)
        menubar.add_cascade(label="Image", menu=image_menu)

        view_menu = tk.Menu(menubar, tearoff=0)
        view_menu.add_command(label="Zoom In", command=lambda: self.zoom_step(1),
                               accelerator="Ctrl++")
        view_menu.add_command(label="Zoom Out", command=lambda: self.zoom_step(-1),
                               accelerator="Ctrl+-")
        view_menu.add_command(label="Fit to Window", command=self.zoom_fit,
                               accelerator="Ctrl+0")
        menubar.add_cascade(label="View", menu=view_menu)

        ocr_menu = tk.Menu(menubar, tearoff=0)
        ocr_menu.add_command(label="Run OCR", command=self.run_ocr, accelerator="Ctrl+R")
        ocr_menu.add_command(label="Force Reprocess Current Image", command=self.force_reprocess)
        menubar.add_cascade(label="OCR", menu=ocr_menu)

        visit_menu = tk.Menu(menubar, tearoff=0)
        visit_menu.add_command(label="Add Blank Field", command=self.add_blank_field)
        visit_menu.add_command(label="Save Visit && Next Image", command=self.save_and_next,
                                accelerator="Ctrl+S")
        visit_menu.add_command(label="Skip Image (No Save)", command=self.next_image,
                                accelerator="Ctrl+Right")
        menubar.add_cascade(label="Visit", menu=visit_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self.show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menubar)
        self.menubar = menubar

        # Keyboard shortcuts matching the accelerators shown above.
        self.bind_all("<Control-o>", lambda e: self.select_folder())
        self.bind_all("<Control-q>", lambda e: self.destroy())
        self.bind_all("<Control-k>", lambda e: self.auto_scan_current())
        self.bind_all("<Control-r>", lambda e: self.run_ocr())
        self.bind_all("<Control-s>", lambda e: self.save_and_next())
        self.bind_all("<Control-Right>", lambda e: self.next_image())
        self.bind_all("<Control-plus>", lambda e: self.zoom_step(1))
        self.bind_all("<Control-KP_Add>", lambda e: self.zoom_step(1))
        self.bind_all("<Control-minus>", lambda e: self.zoom_step(-1))
        self.bind_all("<Control-KP_Subtract>", lambda e: self.zoom_step(-1))
        self.bind_all("<Control-0>", lambda e: self.zoom_fit())

    def show_about(self):
        messagebox.showinfo(
            "About",
            "Diabetes Records Digitizer\n\n"
            "Scan, OCR, and review FNS follow-up forms into structured "
            "per-patient records.\n\n"
            "All actions are available from the menu bar above: File, "
            "Image, View, OCR, and Visit."
        )

    def _build_ui(self):
        top = ttk.Frame(self)
        top.pack(side=tk.TOP, fill=tk.X, padx=6, pady=4)
        self.progress_label = ttk.Label(top, text="No folder selected")
        self.progress_label.pack(side=tk.LEFT, padx=12)
        self.status_label = ttk.Label(top, text="", foreground="#555555")
        self.status_label.pack(side=tk.RIGHT)

        main = ttk.Frame(self)
        main.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=6, pady=4)

        right_outer = ttk.Frame(main, width=540)
        right_outer.pack(side=tk.RIGHT, fill=tk.BOTH)
        right_outer.pack_propagate(False)

        right_vsb = ttk.Scrollbar(right_outer, orient="vertical")
        right_vsb.pack(side=tk.RIGHT, fill=tk.Y)

        self.right_canvas = tk.Canvas(right_outer, borderwidth=0, highlightthickness=0)
        self.right_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.right_canvas.configure(yscrollcommand=right_vsb.set)
        right_vsb.config(command=self.right_canvas.yview)

        right = ttk.Frame(self.right_canvas)
        self.right_canvas_window = self.right_canvas.create_window((0, 0), window=right, anchor="nw")

        def _on_right_configure(e):
            self.right_canvas.configure(scrollregion=self.right_canvas.bbox("all"))
        right.bind("<Configure>", _on_right_configure)

        def _on_right_canvas_configure(e):
            self.right_canvas.itemconfig(self.right_canvas_window, width=e.width)
        self.right_canvas.bind("<Configure>", _on_right_canvas_configure)

        left = ttk.Frame(main)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        controls = ttk.Frame(left)
        controls.pack(side=tk.TOP, fill=tk.X, pady=4)
        ttk.Label(controls, text="Fine angle:").pack(side=tk.LEFT, padx=(0, 2))
        angle_scale = ttk.Scale(controls, from_=-15, to=15, variable=self.angle_var,
                                 orient=tk.HORIZONTAL, length=150,
                                 command=lambda e: self.on_angle_change())
        angle_scale.pack(side=tk.LEFT)

        ttk.Label(controls, text="Zoom:").pack(side=tk.LEFT, padx=(12, 2))
        self.zoom_label = ttk.Label(controls, text="100%", width=5, anchor="center")
        self.zoom_label.pack(side=tk.LEFT)

        ttk.Label(
            controls,
            text="(drag = crop, wheel = zoom, middle-drag = pan — Scan, Rotate, Denoise, "
                 "Contrast, Shadows, Grid, and Zoom are all in the Image/View menus)"
        ).pack(side=tk.LEFT, padx=8)

        canvas_frame = ttk.Frame(left)
        canvas_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        canvas_frame.rowconfigure(0, weight=1)
        canvas_frame.columnconfigure(0, weight=1)

        self.canvas = tk.Canvas(canvas_frame, bg="#222222", cursor="cross")
        self.canvas.grid(row=0, column=0, sticky="nsew")
        vbar = ttk.Scrollbar(canvas_frame, orient="vertical", command=self.canvas.yview)
        vbar.grid(row=0, column=1, sticky="ns")
        hbar = ttk.Scrollbar(canvas_frame, orient="horizontal", command=self.canvas.xview)
        hbar.grid(row=1, column=0, sticky="ew")
        self.canvas.configure(yscrollcommand=vbar.set, xscrollcommand=hbar.set)

        self.canvas.bind("<ButtonPress-1>", self.on_crop_start)
        self.canvas.bind("<B1-Motion>", self.on_crop_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_crop_end)

        self.canvas.bind("<MouseWheel>", self.on_mousewheel_zoom)
        self.canvas.bind("<Button-4>", self.on_mousewheel_zoom)
        self.canvas.bind("<Button-5>", self.on_mousewheel_zoom)

        self.canvas.bind("<ButtonPress-2>", lambda e: self.canvas.scan_mark(e.x, e.y))
        self.canvas.bind("<B2-Motion>", lambda e: self.canvas.scan_dragto(e.x, e.y, gain=1))

        self.bind_all("<MouseWheel>", self.on_global_mousewheel)
        self.bind_all("<Button-4>", self.on_global_mousewheel)
        self.bind_all("<Button-5>", self.on_global_mousewheel)

        patient_box = ttk.LabelFrame(right, text="Patient & visit")
        patient_box.pack(side=tk.TOP, fill=tk.X, padx=4, pady=4)

        ttk.Label(patient_box, text="Patient ID:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.patient_id_var = tk.StringVar()
        self.patient_combo = ttk.Combobox(patient_box, textvariable=self.patient_id_var, values=[])
        self.patient_combo.grid(row=0, column=1, sticky="we", padx=4, pady=2)
        self.patient_combo.bind("<<ComboboxSelected>>", self.on_patient_selected)
        self.patient_combo.bind("<FocusOut>", self.on_patient_selected)
        self.patient_combo.bind("<Return>", self.on_patient_selected)

        ttk.Label(patient_box, text="Patient name:").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.patient_name_var = tk.StringVar()
        ttk.Entry(patient_box, textvariable=self.patient_name_var).grid(row=1, column=1, sticky="we", padx=4, pady=2)

        ttk.Label(patient_box, text="Record (visit):").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        self.record_var = tk.StringVar()
        self.record_combo = ttk.Combobox(patient_box, textvariable=self.record_var,
                                          values=["+ New record"], state="readonly")
        self.record_combo.set("+ New record")
        self.record_combo.grid(row=2, column=1, sticky="we", padx=4, pady=2)
        self.record_combo.bind("<<ComboboxSelected>>", self.on_record_selected)
        ttk.Label(patient_box, text="(pick an existing record to add this image to it)",
                  foreground="#777777").grid(row=3, column=0, columnspan=2, sticky="w", padx=4)

        ttk.Label(patient_box, text="Visit date:").grid(row=4, column=0, sticky="w", padx=4, pady=2)
        self.visit_date_var = tk.StringVar()
        ttk.Entry(patient_box, textvariable=self.visit_date_var).grid(row=4, column=1, sticky="we", padx=4, pady=2)
        ttk.Label(patient_box, text="(YYYY-MM-DD, auto-suggested from OCR)",
                  foreground="#777777").grid(row=5, column=0, columnspan=2, sticky="w", padx=4)
        patient_box.columnconfigure(1, weight=1)

        review_label = ttk.Label(right, text="Detected fields (edit / label / keep):")
        review_label.pack(side=tk.TOP, anchor="w", padx=4, pady=(8, 0))

        review_container = ttk.Frame(right)
        review_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.review_canvas = tk.Canvas(review_container, borderwidth=0, height=400)
        vsb = ttk.Scrollbar(review_container, orient="vertical", command=self.review_canvas.yview)
        self.review_frame = ttk.Frame(self.review_canvas)
        self.review_frame.bind(
            "<Configure>",
            lambda e: self.review_canvas.configure(scrollregion=self.review_canvas.bbox("all"))
        )
        self.review_canvas.create_window((0, 0), window=self.review_frame, anchor="nw")
        self.review_canvas.configure(yscrollcommand=vsb.set)
        self.review_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        row_header = ttk.Frame(self.review_frame)
        row_header.pack(side=tk.TOP, fill=tk.X)
        for text, w in [("Keep", 4), ("Label", 22), ("Value", 12), ("Unit", 12), ("Conf", 5)]:
            ttk.Label(row_header, text=text, width=w, font=("TkDefaultFont", 9, "bold")).pack(side=tk.LEFT, padx=1)

        bottom_right = ttk.Frame(right)
        bottom_right.pack(side=tk.BOTTOM, fill=tk.X, padx=4, pady=6)
        ttk.Label(
            bottom_right,
            text="Use the Visit menu to add a blank field, save this visit, or skip this image.",
            foreground="#777777"
        ).pack(side=tk.LEFT)

    def on_global_mousewheel(self, event):
        try:
            widget = self.winfo_containing(event.x_root, event.y_root)
        except KeyError:
            return
        if not widget:
            return

        if event.num == 4:
            delta = 1
        elif event.num == 5:
            delta = -1
        else:
            delta = int(event.delta / 120)

        widget_str = str(widget)

        if widget_str.startswith(str(self.canvas)):
            cx = event.x_root - self.canvas.winfo_rootx()
            cy = event.y_root - self.canvas.winfo_rooty()
            direction = 1 if delta > 0 else -1
            self._zoom_at(cx, cy, direction)
        elif widget_str.startswith(str(self.review_canvas)):
            self.review_canvas.yview_scroll(-1 * delta, "units")
        elif widget_str.startswith(str(self.right_canvas)) or any(widget_str.startswith(str(child)) for child in self.right_canvas.children.values()):
            self.right_canvas.yview_scroll(-1 * delta, "units")

    def select_folder(self):
        folder = filedialog.askdirectory(title="Select folder with patient photos")
        if not folder:
            return
        self.store = DataStore(folder)
        all_images = []
        for root, _dirs, files in os.walk(folder):
            if os.path.basename(root) == "patients":
                continue
            for fn in files:
                if os.path.splitext(fn)[1].lower() in IMAGE_EXTENSIONS:
                    all_images.append(os.path.join(root, fn))
        all_images.sort()
        self.image_paths = all_images
        self.patient_combo["values"] = self.store.patient_ids()
        self.current_index = -1
        if not all_images:
            messagebox.showinfo("No images found", "No image files were found in that folder (or subfolders).")
            return
        self.next_image()

    def next_image(self):
        self.current_index += 1
        while (self.current_index < len(self.image_paths)
               and self.store.is_processed(self.image_paths[self.current_index])):
            self.current_index += 1
        if self.current_index >= len(self.image_paths):
            self.progress_label.config(text=f"All {len(self.image_paths)} images processed.")
            messagebox.showinfo("Done", "All images in this folder have been processed.")
            return
        self.load_current_image()

    def force_reprocess(self):
        if self.current_index < 0 or self.current_index >= len(self.image_paths):
            return
        path = self.image_paths[self.current_index]
        h = file_hash(path)
        entry = self.store.manifest["processed"].get(h)
        if entry:
            self.pending_preselect = (entry["patient_id"], entry["visit_index"])
        self.store.manifest["processed"].pop(h, None)
        self.load_current_image()

    def load_current_image(self):
        path = self.image_paths[self.current_index]
        try:
            self.pipeline = ImagePipeline(path)
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            self.next_image()
            return
        self.angle_var.set(0.0)
        self.denoise_var.set(False)
        self.contrast_var.set(False)
        self.grid_var.set(False)
        self.shadow_var.set(False)
        self.zoom_factor = 1.0
        self.ocr_lines = []
        self.grouped_lines = []
        self._clear_review_rows()
        self.current_visit_index = None
        self.patient_id_var.set("")
        self.patient_name_var.set("")
        self.visit_date_var.set("")
        self.record_combo["values"] = ["+ New record"]
        self.record_combo.set("+ New record")
        total = len(self.image_paths)
        done = len(self.store.manifest["processed"])
        self.progress_label.config(
            text=f"Image {self.current_index + 1} / {total}  |  {done} already saved  |  {os.path.basename(path)}"
        )
        self.refresh_display()

        if self.pending_preselect:
            pid, visit_index = self.pending_preselect
            self.pending_preselect = None
            self.patient_id_var.set(pid)
            self.patient_name_var.set(self.store.patient_name(pid))
            self.refresh_record_options()
            for combo_idx, v_idx in enumerate(self.record_index_map):
                if v_idx == visit_index:
                    self.record_combo.current(combo_idx + 1)
                    self.on_record_selected()
                    break

    def refresh_record_options(self):
        pid = self.patient_id_var.get().strip()
        options = self.store.record_options(pid) if (self.store and pid) else []
        self.record_index_map = [idx for idx, _label in options]
        labels = ["+ New record"] + [label for _idx, label in options]
        self.record_combo["values"] = labels
        self.record_combo.set("+ New record")
        self.current_visit_index = None

    def rotate(self, delta):
        if not self.pipeline:
            return
        self.pipeline.rotation_90 = (self.pipeline.rotation_90 + delta) % 360
        self.refresh_display()

    def on_angle_change(self):
        if not self.pipeline:
            return
        self.pipeline.fine_angle = self.angle_var.get()
        self.refresh_display()

    def clear_crop(self):
        if not self.pipeline:
            return
        self.pipeline.crop_rect = None
        self.refresh_display()

    def reset_pipeline(self):
        if not self.pipeline:
            return
        self.pipeline.reset()
        self.angle_var.set(0.0)
        self.denoise_var.set(False)
        self.contrast_var.set(False)
        self.grid_var.set(False)
        self.shadow_var.set(False)
        self.refresh_display()

    def auto_scan_current(self):
        if not self.pipeline:
            return
        self.status_label.config(text="Scanning...")
        self.update_idletasks()

        found = self.pipeline.auto_scan()

        # Sync the manual controls with what the scan just did, so the
        # checkboxes/sliders reflect the pipeline's actual state.
        self.angle_var.set(0.0)
        self.contrast_var.set(True)
        self.shadow_var.set(True)
        self.zoom_factor = 1.0

        # Old OCR box overlays were positioned for the pre-scan geometry;
        # clear them so nothing is shown misaligned. Any values already
        # typed into the review panel on the right are left untouched.
        self.ocr_lines = []
        self.grouped_lines = []

        self.refresh_display()  # ensure pipeline flags picked up before fit
        self.zoom_fit()

        if found:
            self.status_label.config(
                text="Scan complete: document detected, perspective-corrected, "
                     "auto-rotated, deskewed, and shadows removed."
            )
        else:
            self.status_label.config(
                text="Scan complete: no clear document edges found, so the full "
                     "frame was used (still auto-rotated, deskewed, and shadow-corrected). "
                     "Use Crop/Rotate manually if needed."
            )

    def refresh_display(self):
        if not self.pipeline:
            return
        self.pipeline.denoise = self.denoise_var.get()
        self.pipeline.contrast = self.contrast_var.get()
        self.pipeline.shadow_removal = self.shadow_var.get()
        img = self.pipeline.render()
        self._show_image(img, boxes=self.ocr_lines, row_boxes=self.grouped_lines)

    def _show_image(self, bgr_img, boxes=None, row_boxes=None, keep_view=False):
        self.canvas.update_idletasks()
        cw = max(self.canvas.winfo_width(), 400)
        ch = max(self.canvas.winfo_height(), 400)
        h, w = bgr_img.shape[:2]
        self.fit_scale = max(min(cw / w, ch / h), 0.02)
        scale = max(self.fit_scale * self.zoom_factor, 0.02)
        self.display_scale = scale

        prev_xview = self.canvas.xview() if keep_view else None
        prev_yview = self.canvas.yview() if keep_view else None

        disp_w, disp_h = max(1, int(w * scale)), max(1, int(h * scale))
        self.disp_w, self.disp_h = disp_w, disp_h
        interp = cv2.INTER_AREA if scale <= 1.0 else cv2.INTER_LINEAR
        resized = cv2.resize(bgr_img, (disp_w, disp_h), interpolation=interp)
        rgb = cv2.cvtColor(resized, cv2.COLOR_BGR2RGB)
        pil_img = Image.fromarray(rgb)
        self.tk_image = ImageTk.PhotoImage(pil_img)
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor="nw", image=self.tk_image, tags="img")
        self.canvas.config(scrollregion=(0, 0, disp_w, disp_h))

        if self.grid_var.get():
            step = 40
            for x in range(0, disp_w, step):
                self.canvas.create_line(x, 0, x, disp_h, fill="#444444", width=1, tags="grid")
            for y in range(0, disp_h, step):
                self.canvas.create_line(0, y, disp_w, y, fill="#444444", width=1, tags="grid")

        if boxes:
            for i, line in enumerate(boxes):
                conf = line["confidence"]
                color = "#2ecc71" if conf >= self.CONF_GOOD else ("#f1c40f" if conf >= self.CONF_OK else "#e74c3c")
                pts = []
                for x, y in line["bbox"]:
                    pts.extend([x * scale, y * scale])
                self.canvas.create_polygon(*pts, outline=color, fill="", width=2, tags=(f"box_{i}",))

        if row_boxes:
            for i, row in enumerate(row_boxes):
                pts = []
                for x, y in row["bbox"]:
                    pts.extend([x * scale, y * scale])
                self.canvas.create_polygon(*pts, outline="#3498db", fill="", width=1,
                                            dash=(4, 2), tags=(f"row_{i}",))

        if prev_xview is not None:
            self.canvas.xview_moveto(prev_xview[0])
            self.canvas.yview_moveto(prev_yview[0])

        self.zoom_label.config(text=f"{int(self.zoom_factor * 100)}%")

    def zoom_step(self, direction):
        if not self.pipeline:
            return
        cx = self.canvas.winfo_width() / 2
        cy = self.canvas.winfo_height() / 2
        self._zoom_at(cx, cy, direction)

    def zoom_fit(self):
        if not self.pipeline:
            return
        self.zoom_factor = 1.0
        self.refresh_display()
        self.canvas.xview_moveto(0)
        self.canvas.yview_moveto(0)

    def on_mousewheel_zoom(self, event):
        pass

    def _zoom_at(self, widget_x, widget_y, direction):
        old_scale = self.display_scale
        canvas_x = self.canvas.canvasx(widget_x)
        canvas_y = self.canvas.canvasy(widget_y)
        img_x = canvas_x / old_scale
        img_y = canvas_y / old_scale

        step = 1.25 if direction > 0 else (1 / 1.25)
        new_zoom = self.zoom_factor * step
        new_zoom = max(0.2, min(new_zoom, 8.0))
        if abs(new_zoom - self.zoom_factor) < 1e-3:
            return
        self.zoom_factor = new_zoom

        self.refresh_display()

        new_scale = self.display_scale
        new_canvas_x = img_x * new_scale
        new_canvas_y = img_y * new_scale
        frac_x = (new_canvas_x - widget_x) / max(self.disp_w, 1)
        frac_y = (new_canvas_y - widget_y) / max(self.disp_h, 1)
        frac_x = min(max(frac_x, 0.0), 1.0)
        frac_y = min(max(frac_y, 0.0), 1.0)
        self.canvas.xview_moveto(frac_x)
        self.canvas.yview_moveto(frac_y)

    def on_crop_start(self, event):
        self.crop_start = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
        if self.crop_rect_id:
            self.canvas.delete(self.crop_rect_id)
            self.crop_rect_id = None

    def on_crop_drag(self, event):
        if not self.crop_start:
            return
        if self.crop_rect_id:
            self.canvas.delete(self.crop_rect_id)
        x0, y0 = self.crop_start
        cx, cy = self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)
        self.crop_rect_id = self.canvas.create_rectangle(x0, y0, cx, cy, outline="#00d4ff", width=2)

    def on_crop_end(self, event):
        if not self.crop_start or not self.pipeline:
            return
        x0, y0 = self.crop_start
        x1, y1 = self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)
        self.crop_start = None
        if abs(x1 - x0) < 8 or abs(y1 - y0) < 8:
            return
        scale = self.display_scale or 1.0
        rect = (int(x0 / scale), int(y0 / scale), int(x1 / scale), int(y1 / scale))
        if self.pipeline.crop_rect:
            ox1, oy1, _ox2, _oy2 = self.pipeline.crop_rect
            rect = (rect[0] + ox1, rect[1] + oy1, rect[2] + ox1, rect[3] + oy1)
        self.pipeline.crop_rect = rect
        self.refresh_display()

    def run_ocr(self):
        if not self.pipeline:
            return
        self.status_label.config(text="Running OCR... (first run downloads models, please wait)")
        self.update_idletasks()
        img = self.pipeline.render()

        def worker():
            try:
                lines = self.ocr_engine.run(img)
                self.ocr_queue.put(("ok", lines))
            except Exception as e:
                self.ocr_queue.put(("error", str(e)))

        threading.Thread(target=worker, daemon=True).start()
        self.after(200, self._poll_ocr_queue)

    def _poll_ocr_queue(self):
        try:
            status, payload = self.ocr_queue.get_nowait()
        except queue.Empty:
            self.after(200, self._poll_ocr_queue)
            return
        if status == "error":
            self.status_label.config(text="OCR failed")
            messagebox.showerror("OCR error", payload)
            return
        self.ocr_lines = payload
        self.grouped_lines = group_ocr_lines(payload)
        self.status_label.config(
            text=f"OCR found {len(payload)} boxes, reconstructed into {len(self.grouped_lines)} line(s)"
        )
        self._show_image(self.pipeline.render(), boxes=self.ocr_lines, row_boxes=self.grouped_lines)
        self._clear_review_rows(origin="ocr")
        for row_line in self.grouped_lines:
            raw_label, value, unit = parse_label_value_unit(row_line["text"])
            extra_aliases = self.store.aliases if self.store else None
            label = match_known_label(raw_label, extra_aliases=extra_aliases)
            self._add_review_row(
                text=value or row_line["text"], confidence=row_line["confidence"],
                origin="ocr", label=label, unit=unit, raw_label=raw_label,
            )
        self._suggest_date()

    def _suggest_date(self):
        if self.visit_date_var.get():
            return
        for line in self.grouped_lines:
            d = try_parse_date(line["text"])
            if d:
                self.visit_date_var.set(d)
                return
        for line in self.ocr_lines:
            d = try_parse_date(line["text"])
            if d:
                self.visit_date_var.set(d)
                break

    def _clear_review_rows(self, origin=None):
        keep_rows = []
        for row in self.field_rows:
            if origin is None or row["origin"] == origin:
                row["frame"].destroy()
            else:
                keep_rows.append(row)
        self.field_rows = keep_rows

    def add_blank_field(self):
        self._add_review_row(text="", confidence=None, origin="manual")

    def _add_review_row(self, text, confidence, origin="ocr", label="", unit="", raw_label=""):
        frame = ttk.Frame(self.review_frame)
        frame.pack(side=tk.TOP, fill=tk.X, pady=1)

        keep_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frame, variable=keep_var, width=4).pack(side=tk.LEFT, padx=1)

        label_var = tk.StringVar(value=label)
        label_combo = ttk.Combobox(
            frame, 
            textvariable=label_var, 
            values=LABEL_SUGGESTIONS, 
            height=8, 
            width=20
        )
        label_combo.pack(side=tk.LEFT, padx=1)

        value_var = tk.StringVar(value=text)
        value_entry = ttk.Entry(frame, textvariable=value_var, width=12)
        value_entry.pack(side=tk.LEFT, padx=1)

        unit_var = tk.StringVar(value=unit)
        unit_entry = ttk.Entry(frame, textvariable=unit_var, width=12)
        unit_entry.pack(side=tk.LEFT, padx=1)

        conf_text = f"{confidence:.2f}" if confidence is not None else "--"
        conf_color = "#2ecc71" if (confidence or 0) >= self.CONF_GOOD else (
            "#f1c40f" if (confidence or 0) >= self.CONF_OK else "#e74c3c")
        conf_label = ttk.Label(frame, text=conf_text, width=5, foreground=conf_color)
        conf_label.pack(side=tk.LEFT, padx=1)

        remove_btn = ttk.Button(frame, text="x", width=2,
                                 command=lambda: self._remove_row(frame))
        remove_btn.pack(side=tk.LEFT, padx=1)

        tag_text = "(existing)" if origin == "existing" else ""
        warn_label = ttk.Label(frame, text=tag_text, foreground="#7f8c8d")
        warn_label.pack(side=tk.LEFT, padx=4)

        def on_change(*_args):
            warning = validate_value(label_var.get(), value_var.get())
            warn_label.config(text=warning or tag_text, foreground="#e67e22" if warning else "#7f8c8d")

        label_var.trace_add("write", on_change)
        value_var.trace_add("write", on_change)

        self.field_rows.append({
            "frame": frame, "keep": keep_var, "label": label_var,
            "value": value_var, "unit": unit_var, "confidence": confidence,
            "origin": origin, "raw_label": raw_label,
        })

    def _remove_row(self, frame):
        self.field_rows = [r for r in self.field_rows if r["frame"] is not frame]
        frame.destroy()

    def on_patient_selected(self, _event=None):
        pid = self.patient_id_var.get().strip()
        if pid and self.store:
            self.patient_name_var.set(self.store.patient_name(pid))
        self.refresh_record_options()

    def on_record_selected(self, _event=None):
        choice = self.record_combo.current()
        if choice <= 0:
            self.current_visit_index = None
            self._clear_review_rows(origin="existing")
            return
        visit_index = self.record_index_map[choice - 1]
        pid = self.patient_id_var.get().strip()
        visit = self.store.get_visit(pid, visit_index)
        if not visit:
            return
        self.current_visit_index = visit_index
        self.visit_date_var.set(visit.get("visit_date") or "")
        self._clear_review_rows(origin="existing")
        existing_fields = visit.get("fields", [])
        for field in existing_fields:
            self._add_review_row(
                text=field.get("value", ""),
                confidence=field.get("confidence"),
                origin="existing",
                label=field.get("label", ""),
                unit=field.get("unit", ""),
            )
        n_images = len(visit.get("images", []))
        self.status_label.config(
            text=f"Loaded record with {len(existing_fields)} existing field(s) from {n_images} image(s). "
                 "This image will be added to it."
        )

    def save_and_next(self):
        if not self.store or not self.pipeline:
            return
        pid = self.patient_id_var.get().strip()
        if not pid:
            messagebox.showwarning("Missing patient", "Please enter or select a Patient ID before saving.")
            return
        name = self.patient_name_var.get().strip()
        visit_date = self.visit_date_var.get().strip()

        fields = []
        learned_count = 0
        for row in self.field_rows:
            if not row["keep"].get():
                continue
            label = row["label"].get().strip()
            value = row["value"].get().strip()
            unit = row["unit"].get().strip()
            if not label and not value:
                continue

            # If this came from OCR and the person corrected the label away
            # from what was auto-detected/auto-filled, and the corrected
            # label matches one of our known canonical fields, remember
            # that abbreviation for next time (this lab's own shorthand).
            raw_label = row.get("raw_label", "")
            if raw_label and normalize_label(raw_label) != normalize_label(label):
                canonical = _LABEL_LOOKUP_NORM.get(normalize_label(label))
                if canonical and self.store.learn_alias(raw_label, canonical):
                    learned_count += 1

            warning = validate_value(label, value)
            spec = KNOWN_FIELDS.get(normalize_label(label))
            fields.append({
                "label": label,
                "value": value,
                "confidence": row["confidence"],
                "unit": unit or (spec["unit"] if spec else ""),
                "flag": warning or "",
            })

        self.store.save_or_update_visit(
            pid, name, self.current_visit_index,
            self.image_paths[self.current_index], visit_date, fields
        )
        self.store.save()
        self.patient_combo["values"] = self.store.patient_ids()
        action = "Updated" if self.current_visit_index is not None else "Saved new"
        msg = f"{action} record for patient {pid} ({len(fields)} fields)"
        if learned_count:
            msg += f" - learned {learned_count} new label synonym(s)"
        self.status_label.config(text=msg)
        self.next_image()

    # ------------------------------------------------------- synonyms UI --
    def open_alias_manager(self):
        if not self.store:
            messagebox.showinfo("No folder open", "Select a folder first - synonyms are saved per folder.")
            return

        win = tk.Toplevel(self)
        win.title("Manage Label Synonyms")
        win.geometry("640x480")

        ttk.Label(
            win,
            text="Different labs abbreviate the same test differently.\n"
                 "Built-in synonyms are always available. Entries you add or that\n"
                 "get learned automatically (by correcting a label) are saved to\n"
                 "label_aliases.json in this folder and used for future auto-fill.",
            justify="left"
        ).pack(side=tk.TOP, anchor="w", padx=10, pady=8)

        columns = ("alias", "canonical", "source")
        tree = ttk.Treeview(win, columns=columns, show="headings", height=15)
        tree.heading("alias", text="Alias (as typed/seen)")
        tree.heading("canonical", text="Maps to")
        tree.heading("source", text="Source")
        tree.column("alias", width=200)
        tree.column("canonical", width=260)
        tree.column("source", width=100, anchor="center")
        tree.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=4)

        def refresh_tree():
            tree.delete(*tree.get_children())
            for canonical, synonyms in sorted(BUILTIN_ALIASES.items()):
                for syn in synonyms:
                    tree.insert("", tk.END, iid=f"builtin::{normalize_label(syn)}",
                                values=(syn, canonical, "built-in"))
            for norm, entry in sorted(self.store.aliases.items(), key=lambda kv: kv[1]["example"]):
                tree.insert("", tk.END, iid=f"learned::{norm}",
                            values=(entry.get("example", norm), entry.get("canonical", ""), "learned/added"))

        refresh_tree()

        add_frame = ttk.LabelFrame(win, text="Add a synonym")
        add_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        ttk.Label(add_frame, text="Alias:").grid(row=0, column=0, padx=4, pady=4, sticky="w")
        alias_var = tk.StringVar()
        ttk.Entry(add_frame, textvariable=alias_var, width=20).grid(row=0, column=1, padx=4, pady=4)

        ttk.Label(add_frame, text="Maps to:").grid(row=0, column=2, padx=4, pady=4, sticky="w")
        canonical_var = tk.StringVar()
        canonical_combo = ttk.Combobox(add_frame, textvariable=canonical_var,
                                        values=LABEL_SUGGESTIONS, width=28)
        canonical_combo.grid(row=0, column=3, padx=4, pady=4)

        def add_alias():
            alias = alias_var.get().strip()
            canonical = canonical_var.get().strip()
            if not alias or not canonical:
                messagebox.showwarning("Missing info", "Enter both an alias and what it maps to.", parent=win)
                return
            if canonical not in LABEL_SUGGESTIONS:
                if not messagebox.askyesno(
                    "Not a known field",
                    f'"{canonical}" isn\'t one of the existing label suggestions.\n'
                    "Add it as a synonym target anyway?", parent=win
                ):
                    return
            self.store.learn_alias(alias, canonical)
            self.store.save_aliases()
            alias_var.set("")
            canonical_var.set("")
            refresh_tree()

        ttk.Button(add_frame, text="Add", command=add_alias).grid(row=0, column=4, padx=8, pady=4)

        def remove_selected():
            sel = tree.selection()
            if not sel:
                return
            changed = False
            for iid in sel:
                if iid.startswith("learned::"):
                    self.store.remove_alias(iid.split("::", 1)[1])
                    changed = True
            if changed:
                self.store.save_aliases()
                refresh_tree()

        bottom = ttk.Frame(win)
        bottom.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=8)
        ttk.Button(bottom, text="Remove Selected (learned/added only)",
                   command=remove_selected).pack(side=tk.LEFT)
        ttk.Button(bottom, text="Close", command=win.destroy).pack(side=tk.RIGHT)


def main():
    if PaddleOCR is None:
        print("WARNING: paddleocr is not installed. The app will open, but "
              "'Run OCR' will fail until you run:\n"
              "  pip install paddleocr paddlepaddle")
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()