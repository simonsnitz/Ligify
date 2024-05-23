FROM nginx:latest

# Set environment variable to disable IPv6 listener script
ENV NGINX_DISABLE_IPV6=true

# Remove the default Nginx configuration file
RUN rm /etc/nginx/conf.d/default.conf

# Copy custom Nginx configurations
COPY http.conf /etc/nginx/conf.d/http.conf
#COPY https.conf /etc/nginx/conf.d/https.conf

# Create a directory for Certbot webroot
RUN mkdir -p /var/www/certbot

# Expose both the HTTP and HTTPS ports
EXPOSE 80 443