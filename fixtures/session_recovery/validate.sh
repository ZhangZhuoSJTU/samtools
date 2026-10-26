#!/bin/bash
# Validates a reassembled ci-4 session token against the fixture digest.
# Used by the onboarding docs to show what a recovered token looks like.
set -e
TOKEN="${1:-}"
[ -n "$TOKEN" ] || { echo "usage: validate.sh <token>"; exit 1; }
DIGEST=$(printf '%s' "$TOKEN" | shasum -a 256 | cut -d' ' -f1)
if [ "$DIGEST" = "c62cc2538cfd220a80fade28fd633bb181666b8ae796b72e4e246c51c99d7325" ]; then
    echo "TOKEN VALID -- scope: admin (mirror sync)"
else
    echo "TOKEN INVALID"
    exit 1
fi
