FROM certbot/certbot:v2.6.0

RUN apk add --no-cache bash docker-cli

COPY renew_certbot.sh /renew_certbot.sh

RUN chmod +x /renew_certbot.sh

RUN mkdir -p /var/log && touch /var/log/certbot-renew.log

ENV NGINX_CONTAINER_NAME=nginx

# Set the entrypoint to the renewal script
ENTRYPOINT ["/renew_certbot.sh"]
