#!/usr/bin/env bash
# Install a user LaunchAgent so the ABC run survives Cursor/agent shell teardown.
# Usage: ./install_abc_launchagent.sh [load|unload|status]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
LABEL="com.zarncke.brain-to-values.abc-extended"
PLIST="$HOME/Library/LaunchAgents/${LABEL}.plist"
LOG="/tmp/abc-extended-r3.log"
ACTION="${1:-load}"

mkdir -p "$HOME/Library/LaunchAgents"
cat >"$PLIST" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key>
  <string>${LABEL}</string>
  <key>WorkingDirectory</key>
  <string>${ROOT}</string>
  <key>ProgramArguments</key>
  <array>
    <string>/bin/bash</string>
    <string>-lc</string>
    <string>set -a; source .env; set +a; exec caffeinate -i bash papers/on-hanson-values/run_abc_extended.sh</string>
  </array>
  <key>StandardOutPath</key>
  <string>${LOG}</string>
  <key>StandardErrorPath</key>
  <string>${LOG}</string>
  <key>RunAtLoad</key>
  <false/>
  <key>KeepAlive</key>
  <false/>
</dict>
</plist>
EOF
echo "Wrote $PLIST"

GUIUID="$(id -u)"
case "$ACTION" in
  load)
    launchctl bootout "gui/${GUIUID}" "$PLIST" 2>/dev/null || true
    launchctl bootstrap "gui/${GUIUID}" "$PLIST"
    launchctl kickstart -k "gui/${GUIUID}/${LABEL}"
    echo "Started ${LABEL} (log: ${LOG})"
    ;;
  unload)
    launchctl bootout "gui/${GUIUID}" "$PLIST" 2>/dev/null || true
    echo "Unloaded ${LABEL}"
    ;;
  status)
    launchctl print "gui/${GUIUID}/${LABEL}" 2>/dev/null || echo "not loaded"
    ;;
  *)
    echo "Usage: $0 [load|unload|status]" >&2
    exit 1
    ;;
esac
