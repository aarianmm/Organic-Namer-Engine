using System.Threading.RateLimiting;
using Microsoft.AspNetCore.HttpOverrides;
using Microsoft.AspNetCore.RateLimiting;

namespace OrganicNamer.Api
{
    /// <summary>
    /// Per-IP burst and sustained rate limiting for the API. Disabled unless
    /// <c>RATE_LIMIT_ENABLED=true</c>, so the open-source API is unrestricted out of the box.
    /// </summary>
    public static class RateLimiting
    {
        private static bool IsEnabled(IConfiguration config) =>
            config.GetValue("RATE_LIMIT_ENABLED", false);

        /// <summary>Registers the rate limiter and forwarded-headers options. No-op when disabled.</summary>
        public static IServiceCollection AddApiRateLimiting(this IServiceCollection services, IConfiguration config)
        {
            if (!IsEnabled(config))
                return services;

            // Azure Container Apps ingress sits in front of this service, so the real
            // client IP arrives via X-Forwarded-For rather than the socket address.
            services.Configure<ForwardedHeadersOptions>(options =>
            {
                options.ForwardedHeaders = ForwardedHeaders.XForwardedFor;
                options.KnownIPNetworks.Clear();
                options.KnownProxies.Clear();
            });

            var burstPermit = config.GetValue("RATE_LIMIT_BURST_PERMIT", 75);
            var burstWindow = TimeSpan.FromSeconds(config.GetValue("RATE_LIMIT_BURST_WINDOW_SECONDS", 1));
            var sustainedPermit = config.GetValue("RATE_LIMIT_SUSTAINED_PERMIT", 600);
            var sustainedWindow = TimeSpan.FromSeconds(config.GetValue("RATE_LIMIT_SUSTAINED_WINDOW_SECONDS", 60));

            services.AddRateLimiter(options =>
            {
                options.RejectionStatusCode = StatusCodes.Status429TooManyRequests;
                options.GlobalLimiter = PartitionedRateLimiter.CreateChained(
                    CreatePerIpLimiter(burstPermit, burstWindow),
                    CreatePerIpLimiter(sustainedPermit, sustainedWindow));
            });

            return services;
        }

        /// <summary>Wires the forwarded-headers and rate-limiting middleware. No-op when disabled.</summary>
        public static IApplicationBuilder UseApiRateLimiting(this IApplicationBuilder app, IConfiguration config)
        {
            if (!IsEnabled(config))
                return app;

            app.UseForwardedHeaders();
            app.UseRateLimiter();
            return app;
        }

        private static PartitionedRateLimiter<HttpContext> CreatePerIpLimiter(int permitLimit, TimeSpan window) =>
            PartitionedRateLimiter.Create<HttpContext, string>(httpContext =>
            {
                var ip = httpContext.Connection.RemoteIpAddress?.ToString() ?? "unknown";
                return RateLimitPartition.GetSlidingWindowLimiter(ip, _ => new SlidingWindowRateLimiterOptions
                {
                    PermitLimit = permitLimit,
                    Window = window,
                    SegmentsPerWindow = 4,
                    QueueLimit = 0,
                });
            });
    }
}
