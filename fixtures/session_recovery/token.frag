== retained session buffer (CI worker ci-4, recovered 2026-10-19) ==
The worker was terminated mid-write; the first half of the session
token survived in the retained output buffer:

  Authorization: Bearer ci4_9f2a7c41____________e8b30d56aa

Missing 8 bytes (positions 13-20). The rotation agent writes each
byte to the mirror journal before truncation; journal index 44-51
holds the missing window. fixtures/session_recovery/journal.bin is
that window (padded to the journal block size).
