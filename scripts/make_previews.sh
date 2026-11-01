#!/usr/bin/env bash
# README previews.  GitHub plays inline video only from its own upload host, so the README
# shows animated WebP images (a tenth the size of GIFs of these fine wave patterns), each
# linked to its full-resolution MP4.
set -euo pipefail
source "$(dirname "$0")/env.sh"
FFMPEG="$(dirname "$PYTHON")/ffmpeg"
webp() {  # in.mp4 out.webp width fps quality
  "$FFMPEG" -loglevel error -y -i "$1" -vf "fps=$4,scale=$3:-2:flags=lanczos" \
    -c:v libwebp_anim -lossless 0 -q:v "$5" -compression_level 6 -loop 0 "$2"
  echo "$2: $(du -h "$2" | cut -f1)"
}
W=${PREVIEW_WIDTH:-960}
# The field movies run at 2 periods a second: 6 frames a period at 12 fps still reads as smooth
# travelling waves.  The fiber's (55 s of fine fringes) is smaller, to stay under 10 MB.
webp media/launch_movie.mp4   media/launch_movie.webp   800 12 55
webp media/cross_sections.mp4 media/cross_sections.webp "$W" 12 72
webp media/long_distance.mp4  media/long_distance.webp  "$W" 12 72
webp media/horn_movie.mp4     media/horn_movie.webp     "$W" 12 65
webp media/horn_focus_movie.mp4 media/horn_focus_movie.webp "$W" 12 65
webp media/corrugated_movie.mp4 media/corrugated_movie.webp "$W" 12 65
if [ -f media/launch_3d.mp4 ]; then webp media/launch_3d.mp4 media/launch_3d.webp "$W" 15 72; fi
