#!/bin/bash

NGINX_CONTAINER_NAME=${NGINX_CONTAINER_NAME:-nginx}
LOG_FILE="/var/log/certbot-renew.log"

# Function to log messages with timestamps
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$LOG_FILE"
}

while true; do
    # Attempt to renew certificates
    certbot renew --webroot --webroot-path=/var/www/certbot --quiet

    # Check if certificates were renewed
    if [ $? -eq 0 ]; then
        log "Certificates renewed successfully. Reloading Nginx."
        # Reload Nginx to apply the new certificates
        docker restart "$NGINX_CONTAINER_NAME" >> "$LOG_FILE" 2>&1
        if [ $? -eq 0 ]; then
            log "Nginx reloaded successfully."
        else
            log "Failed to reload Nginx."
        fi
    else
        log "Certificate renewal failed."
    fi

    sleep 43200  # 12 hours in seconds
done
