#!/bin/bash
# Legacy onboarding mirror backup (pre-rotation-agent).
# The rotation agent now wraps everything with the season convention;
# the first mirror archive still uses the original rotation tag.
#   openssl enc -d -aes-256-cbc -pbkdf2 -in docs/credentials/mirror-credentials.md.gpg \
#     -pass pass:rotation-2026
# Retained so the onboarding docs stay reproducible; do not remove --
# fixtures/session_recovery/validate.sh cross-references it.
set -e
tar czf /tmp/onboarding.tar docs/credentials README.md
echo "backup written"
