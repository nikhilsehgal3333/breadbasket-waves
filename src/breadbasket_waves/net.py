"""Polite HTTP downloads, adapted from the days-of-warning study: one User-Agent, at least
2 seconds between requests to one host, and up to 3 retries after HTTP 429. No personal data goes
into any request. Bodies are returned exactly as the server sent them, so the hashes in the lock
file match the files on the server."""
from __future__ import annotations

import time
import urllib.error
import urllib.request
from urllib.parse import urlparse

USER_AGENT = "breadbasket-waves/0.1 (independent research prototype)"
MIN_GAP = 2.0
_last: dict[str, float] = {}


def get(url: str, timeout: int = 180, tries: int = 3) -> bytes:
    """Return the response body. Raise the last error if all tries fail."""
    host = urlparse(url).netloc
    for attempt in range(tries + 1):
        wait = MIN_GAP - (time.time() - _last.get(host, 0.0))
        if wait > 0:
            time.sleep(wait)
        _last[host] = time.time()
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code == 429 and attempt < tries:
                time.sleep(30)
                continue
            raise
    raise RuntimeError("unreachable")
