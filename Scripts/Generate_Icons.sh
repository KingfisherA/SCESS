#!/bin/sh
# Regenerates the raster icons from Assets/Icon/SCESS Icon.svg (needs Inkscape; the .icns needs macOS)
#
#   Assets/Icon/SCESS Icon.png    512x512, the viewer's splash screen image
#   Assets/Icon/SCESS Icon.icns   macOS application icon, used by the ViewerSCESSApp bundle
#
# Usage: Scripts/Generate_Icons.sh   (from any directory)
set -eu

ICON_DIR="$(cd "$(dirname "$0")/../Assets/Icon" && pwd)"
SVG="$ICON_DIR/SCESS Icon.svg"

if command -v inkscape >/dev/null 2>&1; then
    INKSCAPE=inkscape
elif [ -x /Applications/Inkscape.app/Contents/MacOS/inkscape ]; then
    INKSCAPE=/Applications/Inkscape.app/Contents/MacOS/inkscape
else
    echo "Inkscape not found (macOS: brew install --cask inkscape)" >&2
    exit 1
fi

render() { # <size> <output png>
    "$INKSCAPE" "$SVG" --export-type=png --export-width="$1" --export-height="$1" \
        --export-filename="$2" >/dev/null 2>&1
}

render 512 "$ICON_DIR/SCESS Icon.png"
echo "Wrote $ICON_DIR/SCESS Icon.png"

if command -v iconutil >/dev/null 2>&1; then
    WORK_DIR="$(mktemp -d)"
    trap 'rm -rf "$WORK_DIR"' EXIT
    ICONSET="$WORK_DIR/SCESS.iconset"
    mkdir "$ICONSET"
    for SIZE in 16 32 128 256 512; do
        render "$SIZE" "$ICONSET/icon_${SIZE}x${SIZE}.png"
        render "$((SIZE * 2))" "$ICONSET/icon_${SIZE}x${SIZE}@2x.png"
    done
    iconutil --convert icns --output "$ICON_DIR/SCESS Icon.icns" "$ICONSET"
    echo "Wrote $ICON_DIR/SCESS Icon.icns"
else
    echo "iconutil not found (macOS only); skipped the .icns" >&2
fi
